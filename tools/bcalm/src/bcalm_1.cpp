
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
#ifndef OSX
#include <sys/sysinfo.h> // to determine system memory
#endif


#define CXX11THREADS // let's enable that by default. for now, just comment this line if your gcc version is not high enough
#ifdef CXX11THREADS
 #include <thread>
 #include <atomic>
 #include <../../../thirdparty/concurrentqueue.h>
 #include "../../../thirdparty/ThreadPool.h"
#endif

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
size_t numBucket=1<<minSize;
size_t nb_threads=1;
size_t nb_threads_simulate=1;
bool original_algo = false, use_glueing = true;


// timing-related variables

#ifdef CXX11THREADS
void atomic_double_add(std::atomic<double> &d1, double d2) {
      double current = d1.load();
        while (!d1.compare_exchange_weak(current, current + d2))
                ;
}
typedef std::atomic<double> atomic_double;
#else
#define atomic_double_add(d1,d2) d1 += d2;
typedef double atomic_double;
#endif
atomic_double global_wtime_compactions (0), global_wtime_cdistribution (0), global_wtime_add_nodes (0), global_wtime_create_buckets (0), global_wtime_glue (0), global_wtime_foreach_bucket (0), global_wtime_flush_sb (0), global_wtime_lambda (0), global_wtime_parallel (0), global_wtime_longest_lambda (0), global_wtime_best_sched(0);

bool time_lambdas = true;
std::mutex lambda_timing_mutex;


/********************************************************************************/


bcalm_1::bcalm_1 ()  : Tool ("bcalm_1"){
	getParser()->push_front (new OptionOneParam ("-in", "input file",  true));
	getParser()->push_front (new OptionOneParam ("-out", "output file",  false,"out.fa"));
	getParser()->push_front (new OptionOneParam ("-k", "kmer size",  false,"31"));
	getParser()->push_front (new OptionOneParam ("-m", "minimizer size",  false,"8"));
	getParser()->push_front (new OptionOneParam ("-abundance", "abundance threeshold",  false,"3"));
	getParser()->push_front (new OptionOneParam ("-threads", "number of threads",  false,"2")); // todo: set to max, as in dsk
    getParser()->push_front (new OptionNoParam ("-original", "original BCALM 1 algorithm, without reconciliation",  false));
	getParser()->push_front (new OptionOneParam ("-threads-simulate", "(debug) number of threads to compute scheduling for",  false,"2"));
}


// common to glue and non-glue
void insertInto(vector<BankBinary*>& dests, bool isSuperBucket, size_t index, string seq)
{
	Sequence buffer (Data::ASCII);
	buffer.getData().setRef ((char*)seq.c_str(), seq.size());
	if(dests[index]==NULL){
		BankBinary* bb= new BankBinary((isSuperBucket?"SB":"B")+to_string(index));
		dests[index]=bb;
	}
	dests[index]->insert(buffer);
}

// non-glue
void buckerOrSuper(const string& tmp, size_t min, size_t minimizer,vector<BankBinary*>& superBuckets,vector<BankBinary*>& Buckets){
	size_t prefix(minimizer/numBucket);
	size_t Lprefix(min/numBucket);
	if(Lprefix>prefix){
		insertInto(superBuckets, true, Lprefix, tmp);
	}else{
		insertInto(Buckets, false, min%numBucket, tmp);
	}
}


// non-glue
void goodplace(string& seq, size_t minimizer,BankFasta& out,vector<BankBinary*>& superBuckets,vector<BankBinary*>& Buckets,Model& modelK1){
	Model::Kmer kmmerBegin=modelK1.codeSeed(seq.substr(0,kmerSize-1).c_str(),Data::ASCII);
	size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
	Model::Kmer kmmerEnd=modelK1.codeSeed(seq.substr(seq.size()-kmerSize+1,kmerSize-1).c_str(),Data::ASCII);
	size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));

    if(leftMin<=rightMin){
		if(leftMin>minimizer){
			buckerOrSuper(seq,leftMin,minimizer,superBuckets,Buckets);
		}else{
			if(rightMin>minimizer){
				buckerOrSuper(seq,rightMin,minimizer,superBuckets,Buckets);
			}else{ // leftmin <= minimizer and rightmin <= minimizer, means leftmin=rightmin=minimizer
				Sequence s (Data::ASCII);
				s.getData().setRef ((char*)seq.c_str(), seq.size());
				out.insert(s);
			}
		}
	}else{ // leftmin > rightmin
		if(rightMin>minimizer){ // leftmin > rightmin > minimizer
			buckerOrSuper(seq,rightMin,minimizer,superBuckets,Buckets);
		}else{
			if(leftMin>minimizer){ // leftmin > minimizer >= rightmin
				buckerOrSuper(seq,leftMin,minimizer,superBuckets,Buckets);
			} 
			// else same as above, leftmin <= minimizer and rightmin <= minimizer
		}
	}
}


void put_into_glue(string seq, size_t minimizer, Glue & glue, Model& modelK1) {

	size_t k = kmerSize;
	Model::Kmer kmmerBegin = modelK1.codeSeed(seq.substr(0, k - 1).c_str(), Data::ASCII);
	size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
	Model::Kmer kmmerEnd = modelK1.codeSeed(seq.substr(seq.size() - k + 1, k - 1).c_str(), Data::ASCII);
	size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));

	//size_t i(minimizer/numBucket);
	//size_t j(minimizer%numBucket); // following notation from the ismb2015 paper

	bool lmark = minimizer != leftMin;
	bool rmark = minimizer != rightMin; 
	GlueEntry e(seq, lmark, rmark, kmerSize) ;

	//size_t b_pre_major = leftMin/numBucket;
	//size_t b_suf_major = rightMin/numBucket;

	/*const string keyofin = "CACATGTTCACACACACACACACACACACACACACACACACACACACACAC";
	if ((seq.find(keyofin) != std::string::npos) || (seq.find(reversecomplement(keyofin)) != std::string::npos)) {
		cout << "In put or glue, we got " << tostring(e, e.seq) << endl;
	}
	*/

	glue.insert(e, true);
}


void bcalm_1::execute (){
    static char tableASCII[] = {'A', 'C', 'T', 'G'};
    
    string inputFile(getInput()->getStr("-in"));
    BankFasta out (getInput()->getStr("-out"));
    kmerSize=getInput()->getInt("-k");
    size_t abundance=getInput()->getInt("-abundance");
    minSize=getInput()->getInt("-m");
    nb_threads = getInput()->getInt("-threads");
    nb_threads_simulate = getInput()->getInt("-threads-simulate");
    original_algo = getParser()->saw("-original");
    use_glueing = (!original_algo);
    
    if (nb_threads > nb_threads_simulate)
        nb_threads_simulate = nb_threads;

    Model model(kmerSize, minSize);
    Model modelK1(kmerSize-1, minSize);
    numBucket = 1<<minSize;

    /** check if it's a tiny memory machine, e.g. ec2 micro, if so, limit memory during kmer counting (default is 1G) */
#ifndef OSX
    struct sysinfo info;
    sysinfo(&info);
    int total_ram = (int)(((double)info.totalram*(double)info.mem_unit)/1024/1024);
    printf("Total RAM: %d MB\n",total_ram);
#else
    int total_ram = 128*1024;
#endif
    const char * memory_limit = (total_ram < 1500 ? "-max-memory 500" : "");

    bool is_kmercounted = inputFile.substr(inputFile.size()-2,2) == "h5";

    /** kmer counting **/
    auto start_kc=chrono::system_clock::now();
    {   /**DO NOT REMOVE THOSE BRACKETS**/

        if (!is_kmercounted)
        {
            Graph graph = Graph::create ("-in %s -kmer-size %d  -bloom none -out solidKmers.h5  -abundance-min %d -verbose 1 %s", inputFile.c_str(), kmerSize, abundance, memory_limit);
        }
    }
    auto end_kc=chrono::system_clock::now();
    auto waitedFor_kc=end_kc-start_kc;
    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(1);
    cout<<"Kmer-counting wallclock: "<<chrono::duration_cast<chrono::nanoseconds>(waitedFor_kc).count() / unit <<" secs"<<endl;


    /** We set BankBinary buffer. */
    BankBinary::setBufferSize (1000);

    vector<BankBinary*> superBuckets(numBucket);
    auto start_t=chrono::system_clock::now();
    auto start_part=chrono::system_clock::now();
    size_t maxBucket(0);
    unsigned long nbKmers(0);

    /** partitioning into superbuckets **/
    {
        Storage* storage = StorageFactory(STORAGE_HDF5).load ( is_kmercounted ? inputFile.c_str() : "solidKmers.h5");
        
        LOCAL (storage);
        Group& dskGroup = storage->getGroup("dsk");
        string nbPartitionsStrg = dskGroup.getProperty ("nb_partitions");
        size_t nbPartitions = atol (nbPartitionsStrg.c_str());
        typedef Kmer<SPAN>::Count Count;
        Partition<Count>& partition = dskGroup.getPartition<Count> ("solid", nbPartitions);
        
        /**PUT KMER IN SUPERBUCKET **/
        //Here we read kmers, we want to read sequences TODO (R: more precisely.. we want the input to be sequences of arbitrary length? if so, why? input is always the result of dsk i.e. kmers in h5 format..)
        Iterator<Count>* it = partition.iterator();
        LOCAL (it);
        for (it->first(); !it->isDone(); it->next()){
            nbKmers++;
            Kmer<SPAN>::Type current = it->item().value;
            int minimizer=model.getMinimizerValue(current);

            minimizer /= numBucket;

            insertInto(superBuckets, true, minimizer, model.toString (current));
            
            /** in the parallel version, if the minimizers of the left and right (k-1)-mers
                    are different, then we write to both buckets (actually here, superbuckets) */
            if (use_glueing)
            {
                /*
                 // for some reason that code produces buggy minimizers; why??
                 // I'm keeping it here, because it was a really nasty bug, don't want to make it again. something wrong with getKmer() here?
                Data data (Data::BINARY);
                data.set ((char*) &current, kmerSize);
                size_t leftMin(modelK1.getMinimizerValue(modelK1.getKmer(data,0).value()));
                size_t rightMin(modelK1.getMinimizerValue(modelK1.getKmer(data,1).value()));
                printf("minimizers for kmer %s: %d %d, k-1-mers: %s %s\n",model.toString(current).c_str(), leftMin, rightMin, modelK1.toString(modelK1.getKmer(data,0).value()).c_str(), modelK1.toString(modelK1.getKmer(data,1).value()).c_str());
                */

                string seq = model.toString(current);
                Model::Kmer kmmerBegin=modelK1.codeSeed(seq.substr(0,kmerSize-1).c_str(),Data::ASCII);
                size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
                Model::Kmer kmmerEnd=modelK1.codeSeed(seq.substr(seq.size()-kmerSize+1,kmerSize-1).c_str(),Data::ASCII);
                size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));

                if (leftMin != rightMin)
                {
                    size_t max_minimizer = std::max(leftMin, rightMin);
                    max_minimizer /= numBucket;

                    if (max_minimizer != minimizer) // might be in the same superbucket, in that case, don't insert, we'll take care of it when superbucket->bucket expansion
                        insertInto(superBuckets, true, max_minimizer, model.toString (current));
                }
            }
        }
    }
    auto end_part=chrono::system_clock::now();
    auto waitedFor_part=end_part-start_part;
    cout<<"Partitioned " << nbKmers << " solid kmers into " << numBucket << " major buckets in wall-clock "<<chrono::duration_cast<chrono::nanoseconds>(waitedFor_part).count() / unit <<" secs"<<endl;
    
    Glue glue(kmerSize, out);

    double weighted_best_theoretical_speedup_cumul = 0;
    double weighted_best_theoretical_speedup_sum_times = 0;
    double weighted_best_theoretical_speedup = 0;
    double weighted_actual_theoretical_speedup_cumul = 0;
    double weighted_actual_theoretical_speedup_sum_times = 0;
    double weighted_actual_theoretical_speedup = 0;


    auto start_buckets=chrono::system_clock::now();

    /**FOREACH SUPERBUCKET **/
    for(uint i(0);i<numBucket;++i){
        cout<<'-'<<flush;
        if(superBuckets[i]!=NULL){
            auto start_flush_t=get_wtime();
            superBuckets[i]->flush();
            auto end_flush_t=get_wtime();

            atomic_double_add(global_wtime_flush_sb,diff_wtime(start_flush_t, end_flush_t));

            if(superBuckets[i]->getSize()>0){
                auto start_createbucket_t=get_wtime();
                vector<BankBinary*> Buckets(numBucket);
           
                /* expand a superbucket into buckets */ 
                BankBinary::Iterator itBinary (*superBuckets[i]);
                for (itBinary.first(); !itBinary.isDone(); itBinary.next()){
                    string tmp;
                    size_t leftMin(modelK1.getMinimizerValue(modelK1.getKmer(itBinary->getData(),0).value()));
                    size_t rightMin(modelK1.getMinimizerValue(modelK1.getKmer(itBinary->getData(),itBinary->getData().size()-modelK1.getKmerSize()).value()));
                    for (size_t i=0; i< (itBinary->getDataSize()+3)/4; i++){
                        char b = (itBinary->getData()[i] & 0xFF) ;
                        tmp.push_back(tableASCII[(b>>6)&3]);
                        tmp.push_back(tableASCII[(b>>4)&3]);
                        tmp.push_back(tableASCII[(b>>2)&3]);
                        tmp.push_back(tableASCII[(b>>0)&3]);
                    }
                    tmp.resize(itBinary->getDataSize());
                    
                    if (use_glueing)
                    { // with glue, it's simple: we write to both buckets if leftMin != rightMin
                        if (leftMin / numBucket == i)
                        {
                            insertInto(Buckets, false, leftMin%numBucket, tmp);
                        }
                        if (leftMin != rightMin)
                        {
                            if (rightMin / numBucket == i)
                            {
                                insertInto(Buckets, false, rightMin%numBucket, tmp);
                            }
                        }
                    }
                    else
                    {
                        if(leftMin<=rightMin){
                            if((leftMin/numBucket)==i){
                                insertInto(Buckets, false, leftMin%numBucket, tmp);
                            }else{
                                insertInto(Buckets, false, rightMin%numBucket, tmp);
                            }
                        }else{
                            if((rightMin/numBucket)==i){
                                insertInto(Buckets, false, rightMin%numBucket, tmp);
                            }else{
                                insertInto(Buckets, false, leftMin%numBucket, tmp);
                            }
                        }
                    }
                } // end for itBinary

                auto end_createbucket_t=get_wtime();
                atomic_double_add(global_wtime_create_buckets, diff_wtime(start_createbucket_t, end_createbucket_t));

#ifdef CXX11THREADS
                std::vector<std::thread> threads;
                moodycamel::ConcurrentQueue<std::pair<string, size_t> > glue_queue;
                ThreadPool pool(nb_threads);
#else
                int glue_queue;//dummy
#endif


                // parallel_for doesn't look good, it seems to statically partition the range
                //tbb::parallel_for(4, 0, numBucket, ([&Buckets, &glue_queue, &i, &modelK1, &maxBucket, &superBuckets, &parallel, &out ](long j) {

#if 0
                // TEMPORARY CODE to preload largest bucket from a superbucket (useful to debug parallelism)
                vector<string> test_bucket;
                int j = 0, maxj = 0, maxsize = 0;
                for(;j<numBucket;++j){
                    if(Buckets[j]!=NULL){
                        Buckets[j]->flush();
                        if(Buckets[j]->getSize()>0){
                            if (Buckets[j]->getSize() > maxsize)
                            {
                                maxsize = Buckets[j]->getSize();
                                maxj = j;
                            }
                        }
                    }
                }
                j= maxj;
                BankBinary::Iterator itBinarytest (*Buckets[j]);
                for (itBinarytest.first(); !itBinarytest.isDone(); itBinarytest.next()){
                    string tmp;
                    for (size_t i=0; i< (itBinarytest->getDataSize()+3)/4; i++){
                        char b = (itBinarytest->getData()[i] & 0xFF) ;
                        tmp.push_back(tableASCII[(b>>6)&3]);
                        tmp.push_back(tableASCII[(b>>4)&3]);
                        tmp.push_back(tableASCII[(b>>2)&3]);
                        tmp.push_back(tableASCII[(b>>0)&3]);
                    }
                    tmp=tmp.substr(0,itBinarytest->getDataSize());
                    test_bucket.push_back(tmp);
                }
#endif

                // flush all buckets before spawning threads
                 for(uint j(0);j<numBucket;++j){
                    if(Buckets[j]!=NULL){
                        Buckets[j]->flush();
                    }
                 }
                
                std::vector<double> lambda_timings;
                auto start_foreach_bucket_t=get_wtime();

 
                /**FOREACH BUCKET **/
                for(uint j(0);j<numBucket;++j){
                    if(Buckets[j]!=NULL){
                        if(Buckets[j]->getSize()>0){

/* // code for iterating largest bucket
                                for (auto & test_bucket_str: test_bucket) {
                                    string tmp = test_bucket_str;
                                    string seq = tmp;
                                    Model::Kmer kmmerBegin=modelK1.codeSeed(seq.substr(0,kmerSize-1).c_str(),Data::ASCII);
                                    size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
                                    Model::Kmer kmmerEnd=modelK1.codeSeed(seq.substr(seq.size()-kmerSize+1,kmerSize-1).c_str(),Data::ASCII);
                                    size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));
*/


                            auto lambdaCompact = [&Buckets, &glue_queue, &i, j, &modelK1, &maxBucket, &superBuckets, &out, &glue, &lambda_timings]() {
                                //~ graph1 g(kmerSize);
                                /* add nodes to graph */
                                auto start_nodes_t=get_wtime();
                                size_t actualMinimizer((i<<minSize)+j);
                                graph2 g(kmerSize-1,actualMinimizer,minSize);
                                BankBinary::Iterator itBinary (*Buckets[j]);
                                for (itBinary.first(); !itBinary.isDone(); itBinary.next()){
                                    string tmp;
                                    size_t leftMin(modelK1.getMinimizerValue(modelK1.getKmer(itBinary->getData(),0).value()));
                                    size_t rightMin(modelK1.getMinimizerValue(modelK1.getKmer(itBinary->getData(),itBinary->getData().size()-modelK1.getKmerSize()).value()));
                                    for (size_t i=0; i< (itBinary->getDataSize()+3)/4; i++){
                                        char b = (itBinary->getData()[i] & 0xFF) ;
                                        tmp.push_back(tableASCII[(b>>6)&3]);
                                        tmp.push_back(tableASCII[(b>>4)&3]);
                                        tmp.push_back(tableASCII[(b>>2)&3]);
                                        tmp.push_back(tableASCII[(b>>0)&3]);
                                    }
                                    tmp=tmp.substr(0,itBinary->getDataSize());

                                    g.addleftmin(leftMin);
                                    g.addrightmin(rightMin);
                                    g.addvertex(tmp);

                                }
                                auto end_nodes_t=get_wtime();
                                atomic_double_add(global_wtime_add_nodes, diff_wtime(start_nodes_t, end_nodes_t));

                                /* compact graph*/
                                auto start_dbg=get_wtime();
                                g.debruijn();
                                //~ g.compressh(actualMinimizer);
                                g.compress2();
                                auto end_dbg=get_wtime();
                                atomic_double_add(global_wtime_compactions, diff_wtime(start_dbg, end_dbg));

                                /* distribute nodes (to other buckets, or output, or glue) */
                                auto start_cdistribution_t=get_wtime(); 
                                for(uint32_t i(1);i<g.unitigs.size();++i){ // TODO: determine if i(1) is not a bug, why not i(0)?
                                    if(g.unitigs[i].size()!=0){
                                        if (use_glueing)
                                        {
#if defined CXX11THREADS
                                            glue_queue.enqueue(make_pair<string, size_t>((string)(g.unitigs[i]), (size_t)actualMinimizer));
#else
                                            put_into_glue(g.unitigs[i], actualMinimizer, glue, modelK1);
#endif
                                        }
                                        else /* write to another bucket or output */
                                        {
                                            goodplace(g.unitigs[i],actualMinimizer,out,superBuckets,Buckets,modelK1);

                                        }

                                    }
                                }
                                auto end_cdistribution_t=get_wtime();
                                atomic_double_add(global_wtime_cdistribution, diff_wtime(start_cdistribution_t, end_cdistribution_t));

                                /* what was this code for? */ 
                                //~ size_t node_index(0);
                                //~ for(auto it(g.nodes.begin());it!=g.nodes.end();it++){
                                //~ if(it->size()!=0)
                                //~ {
                                //~ int leftmin = g.leftmins[node_index];
                                //~ int rightmin = g.rightmins[node_index];
                                //~ goodplace(*it,leftmin,rightmin,actualMinimizer,out,numBucket,superBuckets,Buckets);
                                //~ }
                                //~ node_index++;
                                //~ }
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

#ifdef CXX11THREADS
                            /* let's not have one thread per bucket.. we pool them now */
                            /*threads.push_back(
                              std::thread(*/
                            if (nb_threads > 1) 
                                pool.enqueue(lambdaCompact); 
                            else
                                lambdaCompact();
#else
                            lambdaCompact();
#endif
                        } // end if bucket non empty

                    } // end if bucket non null

                } // end for each bucket

#ifdef CXX11THREADS
                //for ( auto& thread : threads ){ thread.join(); }
                pool.join();
#endif

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

                        cout <<"\nIn this superbucket," <<endl;
                        cout <<"                  sum of time spent in lambda's: "<< global_wtime_lambda / 1000000 <<" msecs" <<endl;
                        cout <<"                                 longest lambda: "<< longest_lambda / 1000000 <<" msecs" <<endl;
                        cout <<"         tot time of best scheduling of lambdas: "<< tot_time_best_sched_lambda / 1000000 <<" msecs" <<endl;
                        double best_theoretical_speedup =  global_wtime_lambda  / longest_lambda;
                        double actual_theoretical_speedup =  global_wtime_lambda  / tot_time_best_sched_lambda;
                        cout <<"                       best theoretical speedup: "<<  best_theoretical_speedup << "x" <<endl;
                        if (nb_threads_simulate > 1)
                            cout <<"     best theoretical speedup with "<< nb_threads_simulate << " thread(s): "<<  actual_theoretical_speedup << "x" <<endl;

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

                /* cleaning up */
                {
                    for(uint j(0);j<numBucket;++j){
                        if(Buckets[j]!=NULL){
                            delete(Buckets[j]);
                            remove(("B"+to_string(j)).c_str());
                        }
                    }
                    delete(superBuckets[i]);
                    remove(("SB"+to_string(i)).c_str());
                }

                // do the gluing at the end of each superbucket
                if (use_glueing) {
                    auto start_glue_t=get_wtime();

#ifdef CXX11THREADS
                    // put everything into the glue
                    std::pair<string, size_t> glue_elt;
                    while (glue_queue.try_dequeue(glue_elt))
                    {
                        put_into_glue(glue_elt.first, glue_elt.second, glue, modelK1);
                    }
#endif

					glue.glueStorage.updateMemStats();
					glue.glue();

                    auto end_glue_t=get_wtime();
                    atomic_double_add(global_wtime_glue,  diff_wtime(start_glue_t, end_glue_t));
 
				} // end if parallel
            } // end if superbucket non empty
        } // end if superbucket non null
    } // end foreach superbucket
	
    if (use_glueing) {
		cout << "Final glue:\n" << glue.glueStorage.dump() << "*****\n";
    	glue.glueStorage.printMemStats();
		cout << "Glue class took " << glue.getTime() << "seconds.\n";
	}

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
    cout<<"BCALM total wallclock (incl kmer counting): "<<chrono::duration_cast<chrono::nanoseconds>(end_t-start_t).count() / unit <<" secs"<<endl;
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
