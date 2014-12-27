#include <bcalm_1.hpp>
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



using namespace std;

 const size_t SPAN = KSIZE_1;

    /** Shortcuts. */
    typedef Kmer<SPAN>::Type  Type;
    typedef Kmer<SPAN>::Count Count;
    typedef Kmer<SPAN>::ModelCanonical ModelCanon;
    typedef Kmer<SPAN>::ModelMinimizer <ModelCanon> Model;
    size_t kmerSize=31;
    size_t minSize=8;
    size_t numBucket=1<<minSize;
    size_t threads=2;

/********************************************************************************/


bcalm_1::bcalm_1 ()  : Tool ("bcalm_1"){
    getParser()->push_front (new OptionOneParam ("-in", "input file",  true));
    getParser()->push_front (new OptionOneParam ("-out", "output file",  false,"out.fa"));
    getParser()->push_front (new OptionOneParam ("-k", "kmer size",  false,"31"));
    getParser()->push_front (new OptionOneParam ("-m", "minimizer size",  false,"8"));
    getParser()->push_front (new OptionOneParam ("-abundance", "abundance threeshold",  false,"3"));
    getParser()->push_front (new OptionOneParam ("-threads", "number of threads",  false,"2")); // todo: set to max, as in dsk
}


// todo: write in binary instead of ascii
void insertInto(vector<BankBinary*>& dests, bool isSuperBucket, size_t minimizer, string seq)
{
    Sequence buffer (Data::ASCII);
    buffer.getData().setRef ((char*)seq.c_str(), seq.size());
    if(dests[minimizer]==NULL){
        BankBinary* bb= new BankBinary((isSuperBucket?"SB":"B")+to_string(minimizer));
        dests[minimizer]=bb;
    }
    dests[minimizer]->insert(buffer);
}

void buckerOrSuper(const string& tmp, size_t min, size_t minimizer,vector<BankBinary*>& superBuckets,vector<BankBinary*>& Buckets){
    size_t prefix(minimizer/numBucket);
    size_t Lprefix(min/numBucket);
    if(Lprefix>prefix){
        insertInto(superBuckets, true, Lprefix, tmp);
    }else{
        insertInto(Buckets, false, min%numBucket, tmp);
    }
}


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
            }else{
                Sequence s (Data::ASCII);
                s.getData().setRef ((char*)seq.c_str(), seq.size());
                out.insert(s);
            }
        }
    }else{
        if(rightMin>minimizer){
            buckerOrSuper(seq,rightMin,minimizer,superBuckets,Buckets);
        }else{
            if(leftMin>minimizer){
                buckerOrSuper(seq,leftMin,minimizer,superBuckets,Buckets);
            }else{
                Sequence s (Data::ASCII);
                s.getData().setRef ((char*)seq.c_str(), seq.size());
                out.insert(s);
            }
        }
    }
}


void bcalm_1::execute (){
    static char tableASCII[] = {'A', 'C', 'T', 'G'};
    
    string inputFile(getInput()->getStr("-in"));
    BankFasta out (getInput()->getStr("-out"));
    kmerSize=getInput()->getInt("-k");
    size_t abundance=getInput()->getInt("-abundance");
    minSize=getInput()->getInt("-m");
    threads = getInput()->getInt("-threads");
    
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

    /** kmer counting **/
    auto start_kc=chrono::system_clock::now();
    {
        /**DO NOT REMOVE THOSE BRACKETS**/
        Graph graph = Graph::create ("-in %s -kmer-size %d  -bloom none -out solidKmers.h5  -abundance-min %d -verbose 1 %s", inputFile.c_str(), kmerSize, abundance, memory_limit);
    }
    auto end_kc=chrono::system_clock::now();
    auto waitedFor_kc=end_kc-start_kc;
    cout<<"Kmer-counting wallclock: "<<chrono::duration_cast<chrono::seconds>(waitedFor_kc).count()<<" seconds"<<endl;


    /** We set BankBinary buffer. */
    BankBinary::setBufferSize (1000);

    vector<BankBinary*> superBuckets(numBucket);
    auto start=chrono::system_clock::now();
    auto start_part=chrono::system_clock::now();
    size_t maxBucket(0);
    unsigned long nbKmers(0);
    bool parallel = (threads > 1);  

    //parallel = true;
    // forcing glue code even in single thread for now

    /** partitioning into superbuckets **/
    {
        Storage* storage = StorageFactory(STORAGE_HDF5).load ("solidKmers.h5");
        
        LOCAL (storage);
        Group& dskGroup = storage->getGroup("dsk");
        string nbPartitionsStrg = dskGroup.getProperty ("nb_partitions");
        size_t nbPartitions = atol (nbPartitionsStrg.c_str());
        typedef Kmer<SPAN>::Count Count;
        Partition<Count>& partition = dskGroup.getPartition<Count> ("solid", nbPartitions);
        
        /**PUT KMER IN SUPERBUCKET **/
        //Here we read kmers, we want to read sequences TODO (R: todo what?
        Iterator<Count>* it = partition.iterator();
        LOCAL (it);
        for (it->first(); !it->isDone(); it->next()){
            nbKmers++;
            Kmer<SPAN>::Type current = it->item().value;
            int minimizer=model.getMinimizerValue(current);

            minimizer/=numBucket;

            insertInto(superBuckets, true, minimizer, model.toString (current));
            
            /** in the parallel version, if the minimizers of the left and right (k-1)-mers
                    are different, then we write to both buckets */
            if (parallel)
            {
                Data data (Data::BINARY);
                data.set ((char*) &current, kmerSize);
                size_t leftMin(modelK1.getMinimizerValue(modelK1.getKmer(data,0).value()));
                size_t rightMin(modelK1.getMinimizerValue(modelK1.getKmer(data,1).value()));
                if (leftMin != rightMin)
                {
                    size_t max_minimizer = std::max(leftMin, rightMin);
                    max_minimizer /= numBucket;

                    insertInto(superBuckets, true, max_minimizer, model.toString (current));
                }
            }
        }
    }
    auto end_part=chrono::system_clock::now();
    auto waitedFor_part=end_part-start_part;
    cout<<"Partitioned " << nbKmers << " solid kmers into " << numBucket << " major buckets in wall-clock "<<chrono::duration_cast<chrono::seconds>(waitedFor_part).count()<<" seconds"<<endl;
    
    
    /**FOREACH SUPERBUCKET **/
    for(uint i(0);i<numBucket;++i){
        cout<<'-'<<flush;
        if(superBuckets[i]!=NULL){
            superBuckets[i]->flush();
            if(superBuckets[i]->getSize()>0){
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
                
                /**FOREACH BUCKET **/
                for(uint j(0);j<numBucket;++j){
                    if(Buckets[j]!=NULL){
                        Buckets[j]->flush();
                        if(Buckets[j]->getSize()>0){
                            size_t actualMinimizer((i<<minSize)+j);
                            graph2 g(kmerSize-1,actualMinimizer,minSize);
                            //~ graph1 g(kmerSize);
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
                                
                                g.addvertex(tmp);
                                g.addleftmin(leftMin);
                                g.addrightmin(rightMin);
                            }
                            g.debruijn();
                            //~ g.compressh(actualMinimizer);
                            g.compress2();
                            
                            for(uint32_t i(1);i<g.unitigs.size();++i){
                                if(g.unitigs[i].size()!=0){
                                    goodplace(g.unitigs[i],actualMinimizer,out,superBuckets,Buckets,modelK1);
                                }
                            }
                            
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
                            delete(Buckets[j]);
                            remove(("B"+to_string(j)).c_str());
                        }
                    }
                }
                delete(superBuckets[i]);
                remove(("SB"+to_string(i)).c_str());
            }
        }
    }
    
    auto end=chrono::system_clock::now();
    auto waitedFor=end-start;
    cout<<"BCALM wallclock: "<<chrono::duration_cast<chrono::seconds>(waitedFor).count()<<" seconds"<<endl;
    cout<<"Max bucket : "<<maxBucket<<endl;
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
