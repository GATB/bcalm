#include <bcalm_1.hpp>
#include <ograph.h>
#include <iostream>
#include <memory>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <chrono>



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

/********************************************************************************/


bcalm_1::bcalm_1 ()  : Tool ("bcalm_1"){
    getParser()->push_front (new OptionOneParam ("-in", "input file",  true));
    getParser()->push_front (new OptionOneParam ("-out", "output file",  false,"out.fa"));
    getParser()->push_front (new OptionOneParam ("-nks", "abundance threeshold",  false,"3"));
    getParser()->push_front (new OptionOneParam ("-k", "kmer size",  false,"31"));
    getParser()->push_front (new OptionOneParam ("-m", "minimizer size",  false,"8"));
}

void buckerOrSuper(const string& tmp, size_t min, size_t minimizer,vector<BankBinary*>& superBuckets,vector<BankBinary*>& Buckets){
	Sequence seq (Data::ASCII);
	seq.getData().setRef ((char*)tmp.c_str(), tmp.size());
	size_t prefix(minimizer/numBucket);
	size_t Lprefix(min/numBucket);
	if(Lprefix>prefix){
		if(superBuckets[Lprefix]==NULL){
			BankBinary* bb= new BankBinary("SB"+to_string(Lprefix));
			superBuckets[Lprefix]=bb;
		}
		superBuckets[Lprefix]->insert(seq);
	}else{
		if(Buckets[min%numBucket]==NULL){
			BankBinary* bb= new BankBinary("B"+to_string(min%numBucket));
			Buckets[min%numBucket]=bb;
		}
		Buckets[min%numBucket]->insert(seq);
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
	size_t abundance=getInput()->getInt("-nks");
	minSize=getInput()->getInt("-m");
	
	Model model(kmerSize, minSize);
	Model modelK1(kmerSize-1, minSize);
	numBucket = 1<<minSize;
	
	{
		/**KMER COUNTING DO NOT REMOVE THOSE BRACKETS**/
		Graph graph = Graph::create ("-in %s -kmer-size %d  -bloom none -out solidKmers.h5  -abundance %d -verbose 1", inputFile.c_str(), kmerSize,abundance);
	}
	
	vector<BankBinary*> superBuckets(numBucket);
	auto start=chrono::system_clock::now();
	size_t maxBucket(0);
	{
		Storage* storage = StorageFactory(STORAGE_HDF5).load ("solidKmers.h5");
		
		LOCAL (storage);
		Group& dskGroup = storage->getGroup("dsk");
		string nbPartitionsStrg = dskGroup.getProperty ("nb_partitions");
		size_t nbPartitions = atol (nbPartitionsStrg.c_str());
		typedef Kmer<SPAN>::Count Count;
		Partition<Count>& partition = dskGroup.getPartition<Count> ("solid", nbPartitions);
		
		
		
		/**PUT KMER IN SUPERBUCKET **/
		//Here we read kmers, we want to read sequences TODO
		Iterator<Count>* it = partition.iterator();
		LOCAL (it);
		for (it->first(); !it->isDone(); it->next()){
			Kmer<SPAN>::Type current = it->item().value;
			int minimizer=model.getMinimizerValue(current);
			minimizer/=numBucket;
			
			Sequence seq (Data::ASCII);
			string tmp = model.toString (current);
			
			seq.getData().setRef ((char*)tmp.c_str(), model.getKmerSize());
			
			if(superBuckets[minimizer]==NULL){
				BankBinary* bb= new BankBinary("SB"+to_string(minimizer));
				superBuckets[minimizer]=bb;
			}
			superBuckets[minimizer]->insert(seq);
		}
	}
    
    
    /**FOREACH SUPERBUCKET **/
    for(uint i(0);i<numBucket;++i){
		cout<<'-'<<flush;
		if(superBuckets[i]!=NULL){
			superBuckets[i]->flush();
			if(superBuckets[i]->getSize()>0){
				vector<BankBinary*> Buckets(numBucket);
				
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
					Sequence seq (Data::ASCII);
					seq.getData().setRef ((char*)tmp.c_str(),itBinary->getDataSize());
					if(leftMin<=rightMin){
						if((leftMin/numBucket)==i){
							if(Buckets[leftMin%numBucket]==NULL){
								BankBinary* bb= new BankBinary("B"+to_string(leftMin%numBucket));
								Buckets[leftMin%numBucket]=bb;
							}
							Buckets[leftMin%numBucket]->insert(seq);
						}else{
							if(Buckets[rightMin%numBucket]==NULL){
								BankBinary* bb= new BankBinary("B"+to_string(rightMin%numBucket));
								Buckets[rightMin%numBucket]=bb;
							}
							Buckets[rightMin%numBucket]->insert(seq);
						}
					}else{
						if((rightMin/numBucket)==i){
							if(Buckets[rightMin%numBucket]==NULL){
								BankBinary* bb= new BankBinary("B"+to_string(rightMin%numBucket));
								Buckets[rightMin%numBucket]=bb;
							}
							Buckets[rightMin%numBucket]->insert(seq);
						}else{
							if(Buckets[leftMin%numBucket]==NULL){
								BankBinary* bb= new BankBinary("B"+to_string(leftMin%numBucket));
								Buckets[leftMin%numBucket]=bb;
							}
							Buckets[leftMin%numBucket]->insert(seq);
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
	cout<<"Last for "<<chrono::duration_cast<chrono::seconds>(waitedFor).count()<<" seconds"<<endl;
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
