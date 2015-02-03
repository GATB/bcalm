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
#include <../../../thirdparty/concurrentqueue.h>
#include "../../../thirdparty/ThreadPool.h"

#define get_wtime() chrono::system_clock::now()
#define diff_wtime(x,y) chrono::duration_cast<chrono::nanoseconds>(y - x).count()

using namespace std;

const size_t SPAN = KSIZE_2;

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
atomic_double global_wtime_compactions (0), global_wtime_cdistribution (0), global_wtime_add_nodes (0), global_wtime_create_buckets (0), global_wtime_glue (0), global_wtime_foreach_bucket (0), global_wtime_flush_sb (0), global_wtime_lambda (0), global_wtime_parallel (0), global_wtime_longest_lambda (0), global_wtime_best_sched(0);

bool time_lambdas = true;
std::mutex lambda_timing_mutex, active_minimizers_mutex;


/********************************************************************************/


bcalm_1::bcalm_1 ()  : Tool ("bcalm_1"){
	getParser()->add (new OptionOneParam ("-in", "input file",  true));
	getParser()->add (new OptionOneParam ("-out", "output file",  false,"out.fa"));
	getParser()->add (new OptionOneParam ("-k", "kmer size",  false,"31"));
	getParser()->add (new OptionOneParam ("-m", "minimizer size",  false,"8"));
	getParser()->add (new OptionOneParam ("-abundance", "abundance threeshold",  false,"3"));
	getParser()->add (new OptionOneParam ("-threads-simulate", "(debug) number of threads to compute scheduling for",  false,"1"));
	getParser()->add (new OptionOneParam ("-minimizer-type", "use lexicographical minimizers (0) or frequency based (1)",  false,"1"));
	getParser()->add (new OptionOneParam ("-dsk-memory", "max memory for dsk (MB)", false, "1000"));
}

void bcalm_1::execute (){
    static char tableASCII[] = {'A', 'C', 'T', 'G'};

    string inputFile(getInput()->getStr("-in"));
    BankFasta out (getInput()->getStr("-out"));
    kmerSize=getInput()->getInt("-k");
    size_t abundance=getInput()->getInt("-abundance");
    minSize=getInput()->getInt("-m");
    nb_threads = getInput()->getInt("-nb-cores");
    nb_threads_simulate = getInput()->getInt("-threads-simulate");
    int minimizer_type = getInput()->getInt("-minimizer-type");
    int dsk_memory = getInput()->getInt("-dsk-memory");

    if (nb_threads > nb_threads_simulate)
        nb_threads_simulate = nb_threads;

    /** check if it's a tiny memory machine, e.g. ec2 micro, if so, limit memory during kmer counting (default is 1G) */
#ifndef OSX
    struct sysinfo info;
    sysinfo(&info);
    int total_ram = (int)(((double)info.totalram*(double)info.mem_unit)/1024/1024);
    cout << "\nTotal RAM: " << total_ram << " MB\n";
#else
    int total_ram = 128*1024;
#endif
    if (total_ram < 1500)
        dsk_memory = 500;

    bool is_kmercounted = inputFile.substr(inputFile.size()-2,2) == "h5";

    /** kmer counting **/
    auto start_kc=chrono::system_clock::now();
    {   /**DO NOT REMOVE THOSE BRACKETS**/

        if (!is_kmercounted)
        {
            Graph graph = Graph::create ("-in %s -kmer-size %d -minimizer-size %d  -bloom none -out solidKmers.h5  -abundance-min %d -verbose 1 -minimizer-type %d -max-memory %d", inputFile.c_str(), kmerSize, minSize, abundance, minimizer_type, dsk_memory);
        }
    }
    auto end_kc=chrono::system_clock::now();
    auto waitedFor_kc=end_kc-start_kc;
    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(1);
    cout<<"Kmer-counting wallclock: "<<chrono::duration_cast<chrono::nanoseconds>(waitedFor_kc).count() / unit <<" secs"<<endl;

    /*
     *
     * VARIOUS INIT
     *
     */

    /** We set BankBinary buffer. */
    BankBinary::setBufferSize (1000);

    auto start_t=chrono::system_clock::now();
    auto start_part=chrono::system_clock::now();
    size_t maxBucket(0);
    unsigned long nbKmers(0);

    Storage* storage = StorageFactory(STORAGE_HDF5).load ( is_kmercounted ? inputFile.c_str() : "solidKmers.h5");

    LOCAL (storage);
    Group& dskGroup = storage->getGroup("dsk");
    string nbPartitionsStrg = dskGroup.getProperty ("nb_partitions");
    size_t nbPartitions = atol (nbPartitionsStrg.c_str());
    typedef Kmer<SPAN>::Count Count;
    Partition<Count>& partition = dskGroup.getPartition<Count> ("solid", nbPartitions);
    cout << "DSK created " << nbPartitions << " partitions" << endl;

    /** We retrieve the minimizers distribution from the solid kmers storage. */
    Repartitor repart;
    repart.load (dskGroup);

    u_int64_t rg = ((u_int64_t)1 << (2*minSize));

    /* Retrieve frequency of minimizers;
     * actually only used in minimizerMin and minimizerMax */
    uint32_t *freq_order = NULL;

    if (minimizer_type == 1)
    {
        freq_order = new uint32_t[rg];
        Storage::istream is (dskGroup, "minimFrequency");
        is.read ((char*)freq_order, sizeof(uint32_t) * rg);
    }

    Model model(kmerSize, minSize, Kmer<SPAN>::ComparatorMinimizerFrequency(), freq_order);
    Model modelK1(kmerSize-1, minSize,  Kmer<SPAN>::ComparatorMinimizerFrequency(), freq_order);
    Model modelK2(kmerSize-2, minSize,  Kmer<SPAN>::ComparatorMinimizerFrequency(), freq_order);
    Model modelM(minSize, minSize,  Kmer<SPAN>::ComparatorMinimizerFrequency(), freq_order);

    auto minimizerMin = [&repart, &model] (size_t a, size_t b)
    {
        // should replace with whatever minimizer order function we will use with ModelMinimizer
        return (model.compareIntMinimizers(a,b)) ? a : b;
    };

    auto minimizerMax = [&repart, &model] (size_t a, size_t b)
    {
        return (model.compareIntMinimizers(a,b)) ? b : a;
    };

    /*
     *
     * GLUE
     *
     */

    int nb_glues = 2;
    GlueCommander glue_commander(kmerSize, &out, nb_glues, &model);

    double weighted_best_theoretical_speedup_cumul = 0;
    double weighted_best_theoretical_speedup_sum_times = 0;
    double weighted_best_theoretical_speedup = 0;
    double weighted_actual_theoretical_speedup_cumul = 0;
    double weighted_actual_theoretical_speedup_sum_times = 0;
    double weighted_actual_theoretical_speedup = 0;

    auto start_buckets=chrono::system_clock::now();

    /** We create an iterator for progress information. */
    Iterator<int>* itParts = this->createIterator (
            new Range<int>::Iterator (0,nbPartitions-1), nbPartitions, "Iterating DSK partitions"
            );
    LOCAL (itParts);

    // create many queues in place of Buckets
    std::vector<moodycamel::ConcurrentQueue<std::tuple<string,size_t,size_t> > > bucket_queues;
    bucket_queues.resize(rg);

    /*
     *
     * Iteration of partitions
     *
     */

    std::vector<std::set<size_t>> active_minimizers;
    active_minimizers.resize(nbPartitions);

    /* now our vocabulary is: a "DSK partition" == a "partition" == a "super-bucket"; buckets remainwhat they are in bcalm-original*/
    /* a travelling kmer is one that goes to two buckets from different superbuckets */

    /* main thread is going to read kmers from partitions and insert them into queues  */
    /**FOREACH SUPERBUCKET (= partition) **/
    for (itParts->first (); !itParts->isDone(); itParts->next())
    {
        /** Shortcut. */
        size_t p = itParts->item();

        /** We retrieve an iterator on the Count objects of the pth partition. */
        Iterator<Count>* itKmers = partition[p].iterator();
        LOCAL (itKmers);

        cout << "\nPartition " << p << " has " << partition[p].getNbItems() << " kmers" << endl;

        size_t k = kmerSize;

        auto add_to_bucket_queue = [&active_minimizers, &bucket_queues](size_t minimizer, string seq, size_t leftmin, size_t rightmin, int p)
        {
            bucket_queues[minimizer].enqueue(std::make_tuple(seq,leftmin,rightmin));

            if (active_minimizers[p].find(minimizer) == active_minimizers[p].end())
            {
                active_minimizers_mutex.lock();
                active_minimizers[p].insert(minimizer);
                active_minimizers_mutex.unlock();
            }
        };

        std::atomic<long> nb_left_min_diff_right_min;
        nb_left_min_diff_right_min = 0;

        auto insertIntoQueues = [p, &minimizerMax, &minimizerMin, &add_to_bucket_queue, &bucket_queues, &modelK1, &k, &repart, &nb_left_min_diff_right_min](string seq) {
            Model::Kmer kmmerBegin = modelK1.codeSeed(seq.substr(0, k - 1).c_str(), Data::ASCII);
            size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
            Model::Kmer kmmerEnd = modelK1.codeSeed(seq.substr(seq.size() - k + 1, k - 1).c_str(), Data::ASCII);
            size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));

            if (repart(leftMin) == p)
                add_to_bucket_queue(leftMin, seq, leftMin, rightMin, p);

            if (leftMin != rightMin)
            {
                nb_left_min_diff_right_min ++;

                if (repart(rightMin) == p)
                    add_to_bucket_queue(rightMin, seq, leftMin, rightMin, p);

                // handle traveller kmers
                size_t max_minimizer = minimizerMax(leftMin, rightMin);
                size_t min_minimizer = minimizerMin(leftMin, rightMin);
                if (repart(max_minimizer) != repart(min_minimizer))
                {
                    /* I call that a "traveller kmer" */
                    add_to_bucket_queue(max_minimizer, seq, leftMin, rightMin, repart(max_minimizer));

                    // sanity check
                    if (repart(max_minimizer) < repart(min_minimizer))
                    {                printf("wtf? traveller kmer\n");                exit(1);            }
                }
            }

            // sanity check
            if (repart(leftMin) != p && repart(rightMin) != p)
            {                printf("wtf? repart bucket\n");                exit(1);            }

        };

        auto start_createbucket_t=get_wtime();

        auto iteratePartition = [&insertIntoQueues, &model](Count item) {
            Kmer<SPAN>::Type current = item.value;
            string seq = model.toString(current);
            insertIntoQueues(seq);
        };

        /* expand a superbucket by inserting kmers into queues */
        getDispatcher()->iterate (itKmers, iteratePartition);

        // serial version
        /*for (itKmers->first(); !itKmers->isDone(); itKmers->next()) {
            nbKmers++;
            Kmer<SPAN>::Type current = itKmers->item().value;
            string seq = model.toString(current);
            insertIntoQueues(seq);
        }*/

        cout << "Iterated " << partition[p].getNbItems() << " kmers, among them " << nb_left_min_diff_right_min << " has leftmin!=rightmin" << endl;

        auto end_createbucket_t=get_wtime();
        atomic_double_add(global_wtime_create_buckets, diff_wtime(start_createbucket_t, end_createbucket_t));

        ThreadPool pool(nb_threads - 1);

        std::vector<double> lambda_timings;
        auto start_foreach_bucket_t=get_wtime();

        /**FOREACH BUCKET **/
        for(auto actualMinimizer : active_minimizers[p])
        {
            auto lambdaCompact = [&bucket_queues, actualMinimizer, &glue_commander, &maxBucket, &lambda_timings, &repart, &modelK1]() {
                auto start_nodes_t=get_wtime();

                graph4 g(kmerSize-1,actualMinimizer,minSize);
                //~ //graph1 g(kmerSize);

                /* add nodes to graph */
                std::tuple<string,size_t,size_t> bucket_elt;
                while (bucket_queues[actualMinimizer].try_dequeue(bucket_elt))
                {
                    size_t leftMin(std::get<1>(bucket_elt));
                    size_t rightMin(std::get<2>(bucket_elt));

                    g.addleftmin(leftMin);
                    g.addrightmin(rightMin);
                    g.addvertex(std::get<0>(bucket_elt));
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
                for(uint32_t i(0);i<g.unitigs.size();++i){
                    if(g.unitigs[i].size()!=0){
                        string seq = g.unitigs[i].str();
                        //~ string seq = g.unitigs[i];

                        int k = kmerSize;
                        Model::Kmer kmmerBegin = modelK1.codeSeed(seq.substr(0, k - 1).c_str(), Data::ASCII);
                        size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
                        Model::Kmer kmmerEnd = modelK1.codeSeed(seq.substr(seq.size() - k + 1, k - 1).c_str(), Data::ASCII);
                        size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));

                        bool lmark = actualMinimizer != leftMin;
                        bool rmark = actualMinimizer != rightMin;
                        GlueEntry e(seq, lmark, rmark, kmerSize);
                        glue_commander.insert(e);
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

        glue_commander.cleanup_threaded();
        pool.join();


        if (partition[p].getNbItems() == 0)
            continue; // no stats to print here

        /* compute and print timings */
        auto end_foreach_bucket_t=get_wtime();
        {
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
                glue_commander.queues_size();

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

    // stop the glue thread
    glue_commander.stop();

    // check if buckets are indeed empty
    for (int minimizer = 0; minimizer < rg; minimizer++)
    {
        if  (bucket_queues[minimizer].size_approx() != 0)
        {
            printf("WARNING! bucket %d still has non-processed %d elements\n", minimizer, bucket_queues[minimizer].size_approx() );
            std::tuple<string,size_t,size_t> bucket_elt;
            while (bucket_queues[minimizer].try_dequeue(bucket_elt))
                printf("    %s leftmin %d rightmin %d repartleft %d repartright %d repartmin %d\n", std::get<0>(bucket_elt).c_str(), std::get<1>(bucket_elt), std::get<2>(bucket_elt), repart(std::get<1>(bucket_elt)), repart(std::get<2>(bucket_elt)), repart(minimizer));
        }
    }

    cout << "Final glue:\n";
    glue_commander.dump();
    cout << "*****\n";
    glue_commander.printMemStats();

    /* printing some timing stats */
    auto end_t=chrono::system_clock::now();
    cout<<"Buckets compaction and gluing           : "<<chrono::duration_cast<chrono::nanoseconds>(end_t - start_buckets).count() / unit<<" secs"<<endl;
    cout<<"Within that, \n";
    cout <<"     creating buckets from superbuckets: "<< global_wtime_create_buckets / unit <<" secs"<<endl;
    cout <<"      bucket compaction (except gluing): "<< global_wtime_foreach_bucket / unit <<" secs" <<endl;
    cout <<"     within that, \n";
    cout <<"                       flushing superbuckets: "<< global_wtime_flush_sb / unit <<" secs" <<endl;
    cout <<"                   adding nodes to subgraphs: "<< global_wtime_add_nodes / unit <<" secs" <<endl;
    cout <<"     subgraphs constructions and compactions: "<< global_wtime_compactions / unit <<" secs"<<endl;
    cout <<"              compacted nodes redistribution: "<< global_wtime_cdistribution / unit <<" secs"<<endl;
    cout <<"\nglueing: "<< global_wtime_glue / unit <<" secs"<<endl;
    double sum = global_wtime_glue + global_wtime_cdistribution + global_wtime_compactions + global_wtime_add_nodes + global_wtime_flush_sb + global_wtime_create_buckets;
    cout<<"Sum of the above fine-grained timings: "<< sum / unit <<" secs"<<endl;
    cout<<"Discrepancy between sum of fine-grained timings and total wallclock of buckets compactions step: "<< (chrono::duration_cast<chrono::nanoseconds>(end_t-start_buckets).count() - sum ) / unit <<" secs"<<endl;
    cout<<"BCALM total wallclock (excl kmer counting): "<<chrono::duration_cast<chrono::nanoseconds>(end_t-start_t).count() / unit <<" secs"<<endl;
    cout<<"Max bucket : "<<maxBucket<<endl;
    if (time_lambdas)
    {
        cout<<"\n                 Wallclock time spent in parallel section : "<< global_wtime_parallel / unit << " secs"<<endl;
        cout<<"             Best theoretical speedup in parallel section : "<< weighted_best_theoretical_speedup << "x" <<endl;
        cout<<"Best theoretical speedup in parallel section using " << nb_threads_simulate << " threads : "<< weighted_actual_theoretical_speedup << "x" <<endl;
        cout<<"             Sum of longest bucket compaction for each sb : "<< global_wtime_longest_lambda / unit << " secs"<<endl;
        cout<<"                       Sum of best scheduling for each sb : "<< global_wtime_best_sched / unit << " secs"<<endl;
    }

}




/*
// for some reason that code produces buggy minimizers; why??
// I'm keeping it here, because it was a really nasty bug, don't want to make it again. something wrong with getKmer() here?
Data data (Data::BINARY);
data.set ((char*) &current, kmerSize);
size_t leftMin(modelK1.getMinimizerValue(modelK1.getKmer(data,0).value()));
size_t rightMin(modelK1.getMinimizerValue(modelK1.getKmer(data,1).value()));
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
