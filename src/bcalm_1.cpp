#include <assert.h>
#include <bcalm_1.hpp>
#include <glue.hpp>
#include <ograph.h>
#include <iostream>
#include <memory>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <tuple>
#include "binSeq.h"
#ifndef OSX
#include <sys/sysinfo.h> // to determine system memory
#endif

#include <thread>
#include <atomic>
//#include "../thirdparty/concurrentqueue.h" // not using those
//#include "../thirdparty/lockbasedqueue.h"  // see comments in code for explanations
#include "../thirdparty/lockstdqueue.h" 
//#include "../thirdparty/lockstdvector.h" 

#include "../thirdparty/ThreadPool.h"

#define get_wtime() chrono::system_clock::now()
#define diff_wtime(x,y) chrono::duration_cast<chrono::nanoseconds>(y - x).count()



#define BINSEQ // use binseq in buckets instead of string (slower but uses less memory)
#ifdef BINSEQ
#define BUCKET_STR_TYPE binSeq
#define TO_BUCKET_STR(x) binSeq(x)
#define FROM_BUCKET_STR(x) (x.str())
#else
#define BUCKET_STR_TYPE string
#define TO_BUCKET_STR(x) x
#define FROM_BUCKET_STR(x) x
#endif

using namespace std;

const size_t SPAN = KMER_SPAN(1); // TODO: adapt span using Minia's technique

/** Shortcuts. */
typedef Kmer<SPAN>::Type  Type;
typedef Kmer<SPAN>::Count Count;
typedef Kmer<SPAN>::ModelCanonical ModelCanon;
typedef Kmer<SPAN>::ModelMinimizer <ModelCanon> Model;
size_t kmerSize=31;
size_t minSize=8;
size_t nb_threads=1;
size_t nb_threads_simulate=1;


// timing-related variables

void atomic_double_add(std::atomic<double> &d1, double d2) {
      double current = d1.load();
        while (!d1.compare_exchange_weak(current, current + d2))
                ;
}
typedef std::atomic<double> atomic_double;
#if 0
#define atomic_double_add(d1,d2) d1 += d2;
typedef double atomic_double;
#endif
atomic_double global_wtime_compactions (0), global_wtime_cdistribution (0), global_wtime_add_nodes (0), global_wtime_create_buckets (0), global_wtime_glue_overheads (0), global_wtime_foreach_bucket (0), global_wtime_lambda (0), global_wtime_parallel (0), global_wtime_longest_lambda (0), global_wtime_best_sched(0);

bool time_lambdas = true;
std::mutex lambda_timing_mutex, active_minimizers_mutex, write_to_glue_mutex;

#define time_glue_overhead_start_2() start_glue_t = get_wtime();
#define time_glue_overhead_start() std::chrono::time_point<std::chrono::system_clock> start_glue_t; time_glue_overhead_start_2();
#define time_glue_overhead_stop_2()  end_glue_t=get_wtime(); atomic_double_add(global_wtime_glue_overheads, diff_wtime(start_glue_t, end_glue_t));
#define time_glue_overhead_stop()  std::chrono::time_point<std::chrono::system_clock> end_glue_t; time_glue_overhead_stop_2();

unsigned long memory_usage(string message="")
{
    // using Progress.cpp of gatb-core
    u_int64_t mem = System::info().getMemorySelfUsed() / 1024;
    u_int64_t memMaxProcess = System::info().getMemorySelfMaxUsed() / 1024;
    char tmp[128];
    snprintf (tmp, sizeof(tmp), "  --  memory [current, maximum (maxRSS)]: [%4lu, %4lu] MB ",
            mem, memMaxProcess);
    std::cout << message << " " << tmp << std::endl;
    return mem;
}



/********************************************************************************/


bcalm_1::bcalm_1 ()  : Tool ("bcalm_1"){
	getParser()->push_back (new OptionOneParam ("-in", "input file",  true));
	getParser()->push_back (new OptionOneParam ("-out", "output prefix",  false, "unitigs"));
	getParser()->push_back (new OptionOneParam ("-k", "kmer size",  false,"31"));
	getParser()->push_back (new OptionOneParam ("-m", "minimizer size",  false,"8"));
	getParser()->push_back (new OptionOneParam ("-abundance", "abundance threeshold",  false,"3"));
	getParser()->push_back (new OptionOneParam ("-threads-simulate", "(debug) number of threads to compute scheduling for",  false,"1"));
	getParser()->push_back (new OptionOneParam ("-minimizer-type", "use lexicographical minimizers (0) or frequency based (1)",  false,"1"));
	getParser()->push_back (new OptionOneParam ("-dsk-memory", "max memory for dsk (MB)", false, "1500"));
}

void bcalm_1::execute (){
    static char tableASCII[] = {'A', 'C', 'T', 'G'};

    string inputFile(getInput()->getStr("-in"));
    BankFasta out (getInput()->getStr("-out") + ".fa");
    kmerSize=getInput()->getInt("-k");
    size_t abundance=getInput()->getInt("-abundance");
    minSize=getInput()->getInt("-m");
    nb_threads = getInput()->getInt("-nb-cores");
    nb_threads_simulate = getInput()->getInt("-threads-simulate");
    int minimizer_type = getInput()->getInt("-minimizer-type");
    int dsk_memory = getInput()->getInt("-dsk-memory");

    if (nb_threads > nb_threads_simulate)
        nb_threads_simulate = nb_threads;

    /** check if it's a tiny memory machine, e.g. ec2 micro, if so, limit memory during kmer counting (default is 1.5G) */
#ifndef OSX
    struct sysinfo info;
    sysinfo(&info);
    int total_ram = (int)(((double)info.totalram*(double)info.mem_unit)/1024/1024);
    cout << "\nTotal RAM: " << total_ram << " MB\n";
#else
    int total_ram = 1500; /* don't know how to estimate total ram in osx; maybe GATB knows actually, so todo, use it*/
#endif
    if (total_ram < 1500)
        dsk_memory = 500;

    bool is_kmercounted = inputFile.substr(inputFile.size()-2,2) == "h5";
    string prefix= getInput()->getStr("-out"); 
    if ((!(getInput()->get("-out"))) && is_kmercounted)
    {
        prefix = inputFile.substr(0,inputFile.size()-2);
        cout<< "no -out specified and input is kmercounted: output is  "<< prefix << endl;
    }

    /** kmer counting **/
    bool did_kmercounting = false;
    auto start_kc=chrono::system_clock::now();
    {   /**DO NOT REMOVE THOSE BRACKETS**/

        if (!is_kmercounted)
        {
            Graph graph = Graph::create ("-in %s -kmer-size %d -minimizer-size %d  -bloom none -out %s.h5  -abundance-min %d -verbose 1 -minimizer-type %d -repartition-type 1 -max-memory %d", inputFile.c_str(), kmerSize, minSize, prefix.c_str(), abundance, minimizer_type, dsk_memory);
            did_kmercounting = true;
        }
    }
    auto end_kc=chrono::system_clock::now();
    auto waitedFor_kc=end_kc-start_kc;
    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(1);

    if (did_kmercounting)
    {
        cout<<"Kmer-counting wallclock: "<<chrono::duration_cast<chrono::nanoseconds>(waitedFor_kc).count() / unit <<" secs"<< endl;
        memory_usage("post-kmer counting");
    }

    /*
     *
     * VARIOUS INIT
     *
     */
	initBinSeq();
    /** We set BankBinary buffer. */
    BankBinary::setBufferSize (1000);

    auto start_t=chrono::system_clock::now();
    size_t maxBucket(0);
    unsigned long nbKmers(0);

    Storage* storage = StorageFactory(STORAGE_HDF5).load ( is_kmercounted ? inputFile.c_str() : (prefix + ".h5"));

    LOCAL (storage);
    /** We get the dsk and minimizers hash group in the storage object. */
    Group& dskGroup = storage->getGroup("dsk");
    Group& minimizersGroup = storage->getGroup("minimizers");

    typedef Kmer<SPAN>::Count Count;
    Partition<Count>& partition = dskGroup.getPartition<Count> ("solid");
    size_t nb_partitions = partition.size();
    cout << "DSK created " << nb_partitions << " partitions" << endl;

    /** We retrieve the minimizers distribution from the solid kmers storage. */
    Repartitor repart;
    repart.load (minimizersGroup);

    u_int64_t rg = ((u_int64_t)1 << (2*minSize));

    /* Retrieve frequency of minimizers;
     * actually only used in minimizerMin and minimizerMax */
    uint32_t *freq_order = NULL;

    if (minimizer_type == 1)
    {
        freq_order = new uint32_t[rg];
        Storage::istream is (minimizersGroup, "minimFrequency");
        is.read ((char*)freq_order, sizeof(uint32_t) * rg);
    }

    Model model(kmerSize, minSize, Kmer<SPAN>::ComparatorMinimizerFrequencyOrLex(), freq_order);
    Model modelK1(kmerSize-1, minSize,  Kmer<SPAN>::ComparatorMinimizerFrequencyOrLex(), freq_order);
    Model modelK2(kmerSize-2, minSize,  Kmer<SPAN>::ComparatorMinimizerFrequencyOrLex(), freq_order);
    Model modelM(minSize, minSize,  Kmer<SPAN>::ComparatorMinimizerFrequencyOrLex(), freq_order);

    auto minimizerMin = [&repart, &model] (uint32_t a, uint32_t b)
    {
        return (model.compareIntMinimizers(a,b)) ? a : b;
    };

    auto minimizerMax = [&repart, &model] (uint32_t a, uint32_t b)
    {
        return (model.compareIntMinimizers(a,b)) ? b : a;
    };

    /*
     *
     * GLUE
     *
     */

#ifdef OLD_GLUE
    int nb_glues = 2;
    GlueCommander glue_commander(kmerSize, &out, nb_glues, &model);
#endif

    BankFasta outToGlue (prefix + ".glue");

    double weighted_best_theoretical_speedup_cumul = 0;
    double weighted_best_theoretical_speedup_sum_times = 0;
    double weighted_best_theoretical_speedup = 0;
    double weighted_actual_theoretical_speedup_cumul = 0;
    double weighted_actual_theoretical_speedup_sum_times = 0;
    double weighted_actual_theoretical_speedup = 0;

    auto start_buckets=chrono::system_clock::now();

    /** We create an iterator for progress information. */
    Iterator<int>* it_parts = this->createIterator (
            new Range<int>::Iterator (0,nb_partitions-1), nb_partitions, "Iterating DSK partitions"
            );
    LOCAL (it_parts);


    /*
     *
     * Iteration of partitions
     *
     */

    std::vector<std::set<uint32_t>> active_minimizers;
    active_minimizers.resize(nb_partitions);

    /* now our vocabulary is: a "DSK partition" == a "partition" == a "super-bucket" */
    /* buckets remain what they are in bcalm-original */
    /* a travelling kmer is one that goes to two buckets from different superbuckets */

    // I used to save traveller kmers into bucket_queues, but this would be a memory hog. Let's use files instead. Total volume will be small (a few gigs for human), but that's memory saved
    std::vector<BankFasta*> traveller_kmers_files(nb_partitions);
    std::string traveller_kmers_prefix = prefix + ".doubledKmers.";
    std::mutex *traveller_kmers_save_mutex = new std::mutex[nb_partitions];
    for (unsigned int i = 0; i < nb_partitions; i++)
        traveller_kmers_files[i] = new BankFasta(traveller_kmers_prefix + std::to_string(i));

    auto save_traveller_kmer = [&traveller_kmers_files, &traveller_kmers_save_mutex]
        (uint32_t minimizer, string seq, uint32_t leftmin, uint32_t rightmin, int p) {
            Sequence s (Data::ASCII);
            s.getData().setRef ((char*)seq.c_str(), seq.size());
            traveller_kmers_save_mutex[p].lock();
            traveller_kmers_files[p]->insert(s);
            traveller_kmers_files[p]->flush();
            traveller_kmers_save_mutex[p].unlock();

        };

    /* main thread is going to read kmers from partitions and insert them into queues  */
    /**FOREACH SUPERBUCKET (= partition) **/
    for (it_parts->first (); !it_parts->isDone(); it_parts->next())
    {
        uint32_t p = it_parts->item(); /* partition index */

        /** We retrieve an iterator on the Count objects of the pth partition. */
        Iterator<Count>* it_kmers = partition[p].iterator();
        LOCAL (it_kmers);

        cout << "\nPartition " << p << " has " << partition[p].getNbItems() << " kmers" << endl;

        size_t k = kmerSize;

        // create many queues in place of Buckets
        // (this code used to be outside the partition loop, but I think it's a good idea to reinit the queues after each superbucket(=partition) to avoid queues leaking memory

        // this implementation is supposedly efficient, but:
        // - as fast as the lockbasedqueue below
        // - uses much more memory
        //std::vector<moodycamel::ConcurrentQueue<std::tuple<string,uint32_t,uint32_t> >* > bucket_queues;
        //for (unsigned int i = 0; i < rg; i++)
        //    bucket_queues.push_back(new moodycamel::ConcurrentQueue<std::tuple<string,uint32_t,uint32_t> >);

        // another queue system, very simple, with locks
        // it's fine but uses a linked list, so more memory than I'd like
        //std::vector<LockBasedQueue<std::tuple<string,uint32_t,uint32_t> >* > bucket_queues;
        //for (unsigned int i = 0; i < rg; i++)
        //    bucket_queues.push_back(new LockBasedQueue<std::tuple<string,uint32_t,uint32_t> >);

        // still uses more memory than i'd like
        LockStdQueue<std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t> > bucket_queues[rg]; // TODO: replace string by some sort of binseq class
        //for (unsigned int i = 0; i < rg; i++)
        //    bucket_queues.push_back(new LockStdQueue<std::tuple<string,uint32_t,uint32_t> >);

        //std::vector<LockStdVector<std::tuple<string,uint32_t,uint32_t> >* > bucket_queues;
        //for (unsigned int i = 0; i < rg; i++)
        //    bucket_queues.push_back(new LockStdVector<std::tuple<string,uint32_t,uint32_t> >);


        /* lambda function to add a kmer to a bucket */
        auto add_to_bucket_queue = [&active_minimizers, &bucket_queues](uint32_t minimizer, string seq, uint32_t leftmin, uint32_t rightmin, int p)
        {
            //std::cout << "adding elt to bucket: " << seq << " "<< minimizer<<std::endl;
            bucket_queues[minimizer].enqueue(std::make_tuple(TO_BUCKET_STR(seq),leftmin,rightmin));

            if (active_minimizers[p].find(minimizer) == active_minimizers[p].end())
            {
                active_minimizers_mutex.lock();
                active_minimizers[p].insert(minimizer);
                active_minimizers_mutex.unlock();
            }
        };

        std::atomic<long> nb_left_min_diff_right_min;
        nb_left_min_diff_right_min = 0;
        std::atomic<uint32_t> kmerInGraph;
        kmerInGraph = 0;

        /* lambda function to process a kmer and decide which bucket(s) it should go to */
        auto insertIntoQueues = [p, &minimizerMax, &minimizerMin, &add_to_bucket_queue, 
                    &bucket_queues, &modelK1, &k, &repart, &nb_left_min_diff_right_min,
                    &kmerInGraph, &model, &save_traveller_kmer](Count item) {
            Kmer<SPAN>::Type current = item.value;
            string seq = model.toString(current);

            Model::Kmer kmmerBegin = modelK1.codeSeed(seq.substr(0, k - 1).c_str(), Data::ASCII);
            uint32_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
            Model::Kmer kmmerEnd = modelK1.codeSeed(seq.substr(seq.size() - k + 1, k - 1).c_str(), Data::ASCII);
            uint32_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));

            ++kmerInGraph;

            if (repart(leftMin) == p)
                add_to_bucket_queue(leftMin, seq, leftMin, rightMin, p);

            if (leftMin != rightMin)
            {
                nb_left_min_diff_right_min ++;

                if (repart(rightMin) == p)
                    add_to_bucket_queue(rightMin, seq, leftMin, rightMin, p);

                // handle traveller kmers
                uint32_t max_minimizer = minimizerMax(leftMin, rightMin);
                uint32_t min_minimizer = minimizerMin(leftMin, rightMin);
                if (repart(max_minimizer) != repart(min_minimizer))
                {
                    /* I call that a "traveller kmer" */
                    save_traveller_kmer(max_minimizer, seq, leftMin, rightMin, repart(max_minimizer));
                    //add_to_bucket_queue(max_minimizer, seq, leftMin, rightMin, repart(max_minimizer)); // no longer saved into the queue, but to a file instead

                    // sanity check
                    if (repart(max_minimizer) < repart(min_minimizer))
                    {                printf("wtf? traveller kmer = %s, min_minimizer=%d max_minimizer=%d, repart(min_minimizer)=%d, repart(max_minimizer)=%d\n", seq.c_str(), min_minimizer, max_minimizer, repart(min_minimizer), repart(max_minimizer));                exit(1);            }
                }
            }

            // sanity check
            if (repart(leftMin) != p && repart(rightMin) != p)
            {                printf("wtf? repart bucket\n");                exit(1);            }

        };

        auto start_createbucket_t=get_wtime();

        /* MAIN FIRST LOOP: expand a superbucket by inserting kmers into queues. this creates buckets */
        getDispatcher()->iterate (it_kmers, insertIntoQueues);

        // also add traveller kmers that were saved to a file
        string traveller_kmers_file = traveller_kmers_prefix + std::to_string(p);
        if (System::file().doesExist(traveller_kmers_file)) // for some partitions, there may be no traveller kmers
        {
            BankFasta traveller_kmers_bank (traveller_kmers_file);
            BankFasta::Iterator it (traveller_kmers_bank);
            int nb_traveller_kmers_loaded = 0;
            for (it.first(); !it.isDone(); it.next()) 
            {
                string seq = it->toString();

                // those could be saved in the BankFasta comment eventually
                Model::Kmer kmmerBegin = modelK1.codeSeed(seq.substr(0, k - 1).c_str(), Data::ASCII);
                uint32_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
                Model::Kmer kmmerEnd = modelK1.codeSeed(seq.substr(seq.size() - k + 1, k - 1).c_str(), Data::ASCII);
                uint32_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));

                uint32_t max_minimizer = minimizerMax(leftMin, rightMin);
                add_to_bucket_queue(max_minimizer, seq, leftMin, rightMin, p); 
                nb_traveller_kmers_loaded++;
            }
            std::cout << "Loaded " << nb_traveller_kmers_loaded << " doubled kmers for partition " << p << endl;
            traveller_kmers_bank.finalize();
            System::file().remove (traveller_kmers_file);
        }

        auto end_createbucket_t=get_wtime();
        atomic_double_add(global_wtime_create_buckets, diff_wtime(start_createbucket_t, end_createbucket_t));

        cout << "Iterated " << partition[p].getNbItems() << " kmers, among them " << nb_left_min_diff_right_min << " has leftmin!=rightmin" << endl;

        ThreadPool pool(nb_threads); // FIXME: this used to be "nb_threads-1", but I'm testing without the -1.

        std::vector<double> lambda_timings;
        auto start_foreach_bucket_t=get_wtime();

        /**FOREACH BUCKET **/
        for(auto actualMinimizer : active_minimizers[p])
        {
            auto lambdaCompact = [&bucket_queues, actualMinimizer, 
#ifdef OLD_GLUE
                 &glue_commander, 
#endif
                 &maxBucket, &lambda_timings, &repart, &modelK1, &outToGlue]() {
                auto start_nodes_t=get_wtime();

                // (make sure to change other places labelled "// graph3" and "// graph4" as well)
                //graph4 g(kmerSize-1,actualMinimizer,minSize); // graph4
                graph3 g(kmerSize-1,actualMinimizer,minSize); // graph3
                //~ //graph1 g(kmerSize);

                /* add nodes to graph */
                std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t> bucket_elt;
                while (bucket_queues[actualMinimizer].try_dequeue(bucket_elt))
                {
                    g.addleftmin(std::get<1>(bucket_elt));
                    g.addrightmin(std::get<2>(bucket_elt));
                    g.addvertex(FROM_BUCKET_STR(std::get<0>(bucket_elt)));
                }
                auto end_nodes_t=get_wtime();
                atomic_double_add(global_wtime_add_nodes, diff_wtime(start_nodes_t, end_nodes_t));

                /* compact graph*/
                auto start_dbg=get_wtime();

                g.debruijn();
                g.compress();
                //~ g.compressh(actualMinimizer);

                auto end_dbg=get_wtime();
                atomic_double_add(global_wtime_compactions, diff_wtime(start_dbg, end_dbg));

                /* distribute nodes (to other buckets, or output, or glue) */
                auto start_cdistribution_t=get_wtime();
                 string seq;
                for(uint32_t i(0);i<g.unitigs.size();++i){
                     if(!g.unitigs[i].empty()){ // graph3
    					 seq= g.unitigs[i]; // graph3
					//if(!g.isNumber[i]){ // graph4
						//seq=g.unitigs[i].str(); // graph4


                        int k = kmerSize;
                        Model::Kmer kmmerBegin = modelK1.codeSeed(seq.substr(0, k - 1).c_str(), Data::ASCII);
                        uint32_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
                        Model::Kmer kmmerEnd = modelK1.codeSeed(seq.substr(seq.size() - k + 1, k - 1).c_str(), Data::ASCII);
                        uint32_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));
                        bool lmark = actualMinimizer != leftMin;
                        bool rmark = actualMinimizer != rightMin;

                        bool write_to_disk_instead_of_glueing = true;
                        if (write_to_disk_instead_of_glueing)
                        {
                            if (seq.size() == 0) 
                            {
                                cout<< "compacted unitig smaller than k?! (seq.size() = " << seq.size() << ")" << endl;
                                continue;
                            }
                            Sequence s (Data::ASCII);
                            s.getData().setRef ((char*)seq.c_str(), seq.size());
                            s._comment = string(lmark?"1":"0")+string(rmark?"1":"0"); /** We set the sequence comment. */
                            write_to_glue_mutex.lock();
                            outToGlue.insert(s);
                            outToGlue.flush();
                            write_to_glue_mutex.unlock();
                        }
#ifdef OLD_GLUE
                        else
                        {
                            time_glue_overhead_start();
                            GlueEntry e(seq, lmark, rmark, kmerSize);
                            glue_commander.insert(e);
                            time_glue_overhead_stop();
                        }
#endif
                    }
                }
                auto end_cdistribution_t=get_wtime();
                atomic_double_add(global_wtime_cdistribution, diff_wtime(start_cdistribution_t, end_cdistribution_t));

                size_t size(g.size());
                if(size>maxBucket){
                    maxBucket=size;
                }

                if (time_lambdas)
                {
                    auto time_lambda = diff_wtime(start_nodes_t, end_cdistribution_t);
                    atomic_double_add(global_wtime_lambda, time_lambda);
                    lambda_timing_mutex.lock();
                    lambda_timings.push_back(time_lambda);
                    lambda_timing_mutex.unlock();
                }

            }; // end lambda function

            if (nb_threads > 1)
                pool.enqueue(lambdaCompact);
            else
                lambdaCompact();

        } // end for each bucket

#ifdef OLD_GLUE
        time_glue_overhead_start();
        glue_commander.cleanup_threaded();
        time_glue_overhead_stop();
#endif 
        pool.join();


        if (partition[p].getNbItems() == 0)
            continue; // no stats to print here

        memory_usage("Done with partition " + std::to_string(p));
    
        // check if buckets are indeed empty
        for (unsigned int minimizer = 0; minimizer < rg; minimizer++)
        {
            if  (bucket_queues[minimizer].size_approx() != 0)
            {
                printf("WARNING! bucket %d still has non-processed %d elements\n", minimizer, bucket_queues[minimizer].size_approx() );
                std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t> bucket_elt;
                while (bucket_queues[minimizer].try_dequeue(bucket_elt))
                    printf("    %s leftmin %d rightmin %d repartleft %d repartright %d repartmin %d\n", FROM_BUCKET_STR(std::get<0>(bucket_elt)).c_str(), std::get<1>(bucket_elt), std::get<2>(bucket_elt), repart(std::get<1>(bucket_elt)), repart(std::get<2>(bucket_elt)), repart(minimizer));
            }
        }


        /* compute and print timings */
        {
            auto end_foreach_bucket_t=get_wtime();
            auto wallclock_sb = diff_wtime(start_foreach_bucket_t, end_foreach_bucket_t);
            atomic_double_add(global_wtime_foreach_bucket, wallclock_sb);
            atomic_double_add(global_wtime_parallel, wallclock_sb);

            if (time_lambdas)
            {
                std::sort(lambda_timings.begin(), lambda_timings.end());
                std::reverse(lambda_timings.begin(), lambda_timings.end());
                /* compute a theoretical, i think optimal, scheduling of lambda's using the current number of threads
                */
                double tot_time_best_sched_lambda = 0; // start with the longest lambda
                int t = 0;
                for (auto & lambda_time: lambda_timings) {
                    if ((t++) % nb_threads_simulate == 0)
                        tot_time_best_sched_lambda += lambda_time;
                }

                double longest_lambda = lambda_timings.front();

                cout <<"\nIn this superbucket (containing " << active_minimizers.size() << " active minimizers)," <<endl;
                cout <<"                  sum of time spent in lambda's: "<< global_wtime_lambda / 1000000 <<" msecs" <<endl;
                cout <<"                                 longest lambda: "<< longest_lambda / 1000000 <<" msecs" <<endl;
                cout <<"         tot time of best scheduling of lambdas: "<< tot_time_best_sched_lambda / 1000000 <<" msecs" <<endl;
                double best_theoretical_speedup =  global_wtime_lambda  / longest_lambda;
                double actual_theoretical_speedup =  global_wtime_lambda  / tot_time_best_sched_lambda;
                cout <<"                       best theoretical speedup: "<<  best_theoretical_speedup << "x" <<endl;
                if (nb_threads_simulate > 1)
                    cout <<"     best theoretical speedup with "<< nb_threads_simulate << " thread(s): "<<  actual_theoretical_speedup << "x" <<endl;
                
#ifdef OLD_GLUE
                glue_commander.queues_size();
#endif

                weighted_best_theoretical_speedup_cumul += best_theoretical_speedup * wallclock_sb;
                weighted_best_theoretical_speedup_sum_times                        += wallclock_sb;
                weighted_best_theoretical_speedup = weighted_best_theoretical_speedup_cumul / weighted_best_theoretical_speedup_sum_times ;

                weighted_actual_theoretical_speedup_cumul += actual_theoretical_speedup * wallclock_sb;
                weighted_actual_theoretical_speedup_sum_times                        += wallclock_sb;
                weighted_actual_theoretical_speedup = weighted_actual_theoretical_speedup_cumul / weighted_actual_theoretical_speedup_sum_times ;
                atomic_double_add(global_wtime_longest_lambda, longest_lambda);
                atomic_double_add(global_wtime_best_sched, tot_time_best_sched_lambda);
                global_wtime_lambda = 0;
            }
        }

    } // end iteration superbuckets

    /*
     *
     * Finishing up
     *
     */

#ifdef OLD_GLUE
    // stop the glue thread
    time_glue_overhead_start();
    glue_commander.stop();
    time_glue_overhead_stop();

    cout << "Final glue:\n";
    time_glue_overhead_start_2();
    glue_commander.dump();
    time_glue_overhead_stop_2();
    cout << "*****\n";
    glue_commander.printMemStats();
#endif

    /* printing some timing stats */
    auto end_t=chrono::system_clock::now();
    cout<<"Buckets compaction and gluing           : "<<chrono::duration_cast<chrono::nanoseconds>(end_t - start_buckets).count() / unit<<" secs"<<endl;
    cout<<"Within that, \n";
    cout <<"                                 creating buckets from superbuckets: "<< global_wtime_create_buckets / unit <<" secs"<<endl;
    cout <<"                      bucket compaction (wall-clock during threads): "<< global_wtime_foreach_bucket / unit <<" secs" <<endl;
    cout <<"                                               overhead due to glue: "<< global_wtime_glue_overheads / unit <<" secs"<<endl;
    cout <<"\n                within all bucket compaction threads,\n";
    cout <<"                       adding nodes to subgraphs: "<< global_wtime_add_nodes / unit <<" secs" <<endl;
    cout <<"         subgraphs constructions and compactions: "<< global_wtime_compactions / unit <<" secs"<<endl;
    cout <<"                  compacted nodes redistribution: "<< global_wtime_cdistribution / unit <<" secs"<<endl;
    double sum = global_wtime_glue_overheads + global_wtime_cdistribution + global_wtime_compactions + global_wtime_add_nodes + global_wtime_create_buckets;
    cout<<"Sum of CPU times for bucket compactions: "<< sum / unit <<" secs"<<endl;
    if (nb_threads == 1)
        cout<<"Discrepancy between sum of fine-grained timings and total wallclock of buckets compactions step: "<< (chrono::duration_cast<chrono::nanoseconds>(end_t-start_buckets).count() - sum ) / unit <<" secs"<<endl;
    cout<<"BCALM total wallclock (excl kmer counting): "<<chrono::duration_cast<chrono::nanoseconds>(end_t-start_t).count() / unit <<" secs"<<endl;
    cout<<"Max bucket : "<<maxBucket<<endl;
    if (time_lambdas)
    {
        cout<<"Performance of compaction step:\n"<<endl;
        cout<<"                 Wallclock time spent in parallel section : "<< global_wtime_parallel / unit << " secs"<<endl;
        cout<<"             Best theoretical speedup in parallel section : "<< weighted_best_theoretical_speedup << "x" <<endl;
        cout<<"Best theor. speedup in parallel section using " << nb_threads_simulate << " threads : "<< weighted_actual_theoretical_speedup << "x" <<endl;
        cout<<"             Sum of longest bucket compaction for each sb : "<< global_wtime_longest_lambda / unit << " secs"<<endl;
        cout<<"                       Sum of best scheduling for each sb : "<< global_wtime_best_sched / unit << " secs"<<endl;
    }

}




/*
// for some reason that code produces buggy minimizers; why??
// I'm keeping it here, because it was a really nasty bug, don't want to make it again. something wrong with getKmer() here?
Data data (Data::BINARY);
data.set ((char*) &current, kmerSize);
uint32_t leftMin(modelK1.getMinimizerValue(modelK1.getKmer(data,0).value()));
uint32_t rightMin(modelK1.getMinimizerValue(modelK1.getKmer(data,1).value()));
printf("minimizers for kmer %s: %d %d, k-1-mers: %s %s\n",model.toString(current).c_str(), leftMin, rightMin, modelK1.toString(modelK1.getKmer(data,0).value()).c_str(), modelK1.toString(modelK1.getKmer(data,1).value()).c_str());
*/



 //~ Data data (Data::BINARY);
        //~ data.set ((char*) &current, kmerSize);


            //~ superBuckets[i]->iterate ([&] (const Sequence& s){
            //~ /** We get the kmer corresponding to the current sequence. */
            //~ ModelCanon::Kmer mini = model.codeSeed (s.getDataBuffer(), Data::BINARY);
//~
            //~ cout << "mini=" << mini.value() << "  " << model.toString (mini.value()) << endl;
        //~ });
            //~ cin.get();


                 //~ cout << "[0]=" << modelK1.toString(modelK1.getKmer(itBinary->getData(),0).value())
         //~ << "  mini=" << modelK1.getMmersModel().toString(modelK1.getKmer (itBinary->getData(), 0).minimizer().value() )
         //~ << endl;
//~
//~
    //~ cout << "[N]=" << modelK1.toString(modelK1.getKmer (itBinary->getData(), itBinary->getData().size()-modelK1.getKmerSize()).value())
         //~ << "  mini=" << modelK1.getMmersModel().toString(modelK1.getKmer (itBinary->getData(), itBinary->getData().size()-modelK1.getKmerSize()).minimizer().value() )
         //~ << endl;



        //~ cout << "kmer=" << it->item().value << "  minimizer=" << model.getMinimizerValue(current)
            //~ << "  abundance=" << it->item().abundance << endl;
        //~ cin.get();
        //~ modelK1.iterate (data, [&] (const Model::Kmer& k, size_t idx)
        //~ {
            //~ cout << "-> " << k.value() << " minimizer=" << k.minimizer().value() << "  data.size=" << data.size() << endl;
        //~ });

        //~ Data data (Data::ASCII);
    //~ char* str = "AAAAAAAAAAAAACGTACGATTTTTTTTTTATATAGGATAAAAAAAAAAAACGACTGATCATCGATCATAAAAA";
    //~ data.set (str, strlen(str) );
//~
    //~ cout << "[0]=" << modelK1.toString(modelK1.getKmer (data, 0).value() )
         //~ << "  mini=" << modelK1.getMmersModel().toString(modelK1.getKmer (data, 0).minimizer().value() )
         //~ << endl;
//~
//~
    //~ cout << "[N]=" << modelK1.toString(modelK1.getKmer (data, data.size()-modelK1.getKmerSize()).value())
         //~ << "  mini=" << modelK1.getMmersModel().toString(modelK1.getKmer (data, data.size()-modelK1.getKmerSize()).minimizer().value() )
         //~ << endl;
