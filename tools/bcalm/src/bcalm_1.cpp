
// FIXME: looks buggy on 1000nodechain-k16.fa (compare with bcalm-original, some nodes missing, an extra compaction on the big unitig)
#include <assert.h>
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

//#define SPARSEHASH 
#ifdef SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;
#endif 

//#ifdef _OPENMP
//   #include <omp.h>
//   #define omp_diff_wtime_bcalm(x,y) (y - x) 
//#else
   #define omp_get_thread_num() 0
   #define omp_set_num_threads(x) 0
   #define omp_get_max_threads() 1
   #define omp_get_wtime() chrono::system_clock::now() 
   #define omp_diff_wtime_bcalm(x,y) chrono::duration_cast<chrono::nanoseconds>(y - x).count() 
//#endif


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
size_t threads=1;



template<typename T>
string add_commas(T num) {
	string s, retval;
	ostringstream o;
	o << num;
	s = o.str();
	for (int i = 0; i < s.length(); i++) {
		retval.push_back(s.at(i));
		int j = s.length() - 1 - i;
		if (((j % 3) == 0) && (j != 0)) {
			retval.push_back(',');
		}
	}
	return retval;
}




    // timing-related variables
    double global_wtime_compactions = 0, global_wtime_cdistribution = 0, global_wtime_add_nodes = 0, global_wtime_create_buckets = 0, global_wtime_glue = 0, global_wtime_foreach_bucket = 0, global_wtime_flush_sb = 0, global_wtime_test = 0;
/********************************************************************************/


bcalm_1::bcalm_1 ()  : Tool ("bcalm_1"){
	getParser()->push_front (new OptionOneParam ("-in", "input file",  true));
	getParser()->push_front (new OptionOneParam ("-out", "output file",  false,"out.fa"));
	getParser()->push_front (new OptionOneParam ("-k", "kmer size",  false,"31"));
	getParser()->push_front (new OptionOneParam ("-m", "minimizer size",  false,"8"));
	getParser()->push_front (new OptionOneParam ("-abundance", "abundance threeshold",  false,"3"));
	getParser()->push_front (new OptionOneParam ("-threads", "number of threads",  false,"2")); // todo: set to max, as in dsk
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



string reversecomplement(const string& dna) // code taken from TakeABreak
{
	string revComp= "";

	for (string::const_reverse_iterator it = dna.rbegin(); it != dna.rend(); it++) {
		switch(*it) {
			case 'a' :                revComp += "t";                break;
			case 't' :                revComp += "a";                break;
			case 'c' :                revComp += "g";                break;
			case 'g' :                revComp += "c";                break;
			case 'A' :                revComp += "T";                break;
			case 'T' :                revComp += "A";                break;
			case 'C' :                revComp += "G";                break;
			case 'G' :                revComp += "C";                break;
		}
	}
    return revComp;
}

string int2string(int c)
{
    if (c == 0)
        return "0";
    else if (c == 1)
        return "1";
    else
        return "2";
}

/* nifty function to highlight a substring for printing to terminal */
string debug_highlight(string s, string motif)
{
    size_t pos = s.find(motif);
    if (pos == string::npos)
        pos = s.find(reversecomplement(motif));
    if (pos == string::npos)
    {
        return s;
    }

    return s.substr(0,pos) + "\033[1;31m" + s.substr(pos, motif.size()) + "\033[0m" + s.substr(pos + motif.size(), s.size() - pos - motif.size()); 
}



/*class to support two-bit encoding of entries in the glue table
 * PM: NOT YET ACTIVE
class SeqAndTags {
	public:
		const static int MAX_SIZE = 1000;
		bool raw[MAX_SIZE];
		int rawLen;

		SeqAndTags(bool * _raw, int _rawLen) : rawLen(_rawLen) { 
			for (int i = 0; i < rawLen; i++) {
				raw[i] = _raw[i];
			}
		}

		SeqAndTags(string seq, bool leftmark, bool rightmark) {
			rawLen = seq.length() * 2 + 2;
			for (int i = 0; i < seq.length(); i++) {
				if ((seq[i] == 'A') || (seq[i] == 'C')) {
					raw[2*i] = 0;
				} else {
					raw[2*i] = 1;
				}
				if ((seq[i] == 'A') || (seq[i] == 'G')) {
					raw[2*i + 1] = 0;
				} else {
					raw[2*i + 1] = 1;
				}
			}
			raw[seq.length() * 2] = leftmark;
			raw[seq.length() * 2 + 1] = rightmark;
		}


		char bits2char (bool leftbit, bool rightbit) {
			if (leftbit == 0 && rightbit == 0) return 'A';
			else if (leftbit == 0 && rightbit == 1) return 'C';
			else if (leftbit == 1 && rightbit == 0) return 'G';
			else return 'T';
		}


		string getSeq() {
			string retval;
			for (int  i = 0; i < rawLen - 2; i += 2) {
				retval.push_back(bits2char(raw[i], raw[i+1]));
			}
			return retval;
		}

		string getLeftKmer(int k) {
			return getSeq().substr(0,k);
		}
		string getRightKmer(int k) {
			string seq = getSeq();
			return seq.substr(seq.length() - k, k);
		}
		bool getLeftMark() {
			return raw[rawLen - 2];
		}
		bool getRightMark() {
			return raw[rawLen - 1];
		}
};
*/

class Glue
{
    // optimize it later.. (2 bits sequences)

    public:

    Glue(size_t kmerSize, BankFasta &out) : model(kmerSize), out(out), debug(false)
    {
#ifdef SPARSEHASH
    glueStrings.set_deleted_key("");
#endif
    }

	void now_really_addhash(string key, string value);
	void aux_addhash(string seq, string kmer, bool leftmark, bool rightmark, bool left);
	void insert(string seq, bool leftmark, bool rightmark);
	string seq_only(string seq);
	void remove(string seq);
	void output(string seq);
	void glue();


	void resetMemStats() {
		maxEntries = 0; totEntries = 0; numDataPoints = 0; maxSize = 0; totSize = 0;
	}

	void updateMemStats() {
		maxEntries = std::max(maxEntries, glueStrings.size());
		totEntries += glueStrings.size();

		size_t size = 0;
		for (auto it_gS = glueStrings.begin(); it_gS != glueStrings.end(); it_gS++)
		{
			size += sizeof(it_gS->first) + sizeof(it_gS->second.first) + sizeof(it_gS->second.second);
			size += it_gS->first.size() + it_gS->second.first.size() + it_gS->second.second.size();
		}

		maxSize = std::max(maxSize, size);
		totSize += size;

		numDataPoints++;
	}



	void printMemStats() {
		if (numDataPoints == 0) {
			cout << "Glue: no data points to output memory stats.\n";
		} else {
			cout << "Glue memory stats: max entries: " << add_commas(maxEntries) << " avg entries: " << add_commas(totEntries/numDataPoints) << " max size: " << add_commas(maxSize) 
				<< "b avg size: " << add_commas(totSize / numDataPoints) << "b\n";
		}
	}

	private:

	ModelCanon model;
#ifdef SPARSEHASH
	typedef  sparse_hash_map<string, pair<string, string>> GlueMap; 
#else
	typedef unordered_map<string, pair<string, string>> GlueMap;
#endif
	GlueMap glueStrings;
	BankFasta out;
	bool debug;

	//used for tracking memory usage
	size_t maxEntries = 0;
	size_t totEntries = 0;
	int numDataPoints = 0;
	size_t maxSize = 0;
	size_t totSize = 0;


	
};


    void Glue::now_really_addhash(string key, string value)
    {
          GlueMap::const_iterator got = glueStrings.find (key);
          if ( got == glueStrings.end() )
          {
              glueStrings[key].first = value;
              glueStrings[key].second = "";
          }
          else
          {
              if (glueStrings[key].first == "-1" || glueStrings[key].first == "")
              {
                  glueStrings[key].first = value;
                  glueStrings[key].second = "";
              }
              else
              {
                  if  (glueStrings[key].second != "")
                  {
                      printf("huh? kmer %s (%s,%s) inserting %s\n",key.c_str(), debug_highlight(glueStrings[key].first,key).c_str(), debug_highlight(glueStrings[key].second,key).c_str(), debug_highlight(value,key).c_str());
                      exit(1);
                  }
                  glueStrings[key].second = value;
              }
          }

    }


    void Glue::aux_addhash(string seq, string kmer, bool leftmark, bool rightmark, bool left)
    {
         if (kmer.compare(reversecomplement(kmer)) <= 0)
               now_really_addhash(kmer, seq + int2string(leftmark) + int2string(rightmark));
         else
               now_really_addhash(reversecomplement(kmer), reversecomplement(seq) + int2string(rightmark) + int2string(leftmark));
        return;

        // this was the former strategy
        if (kmer.compare(reversecomplement(kmer)) <= 0)
        {
            if (left)
                glueStrings[kmer].second = seq + int2string(leftmark) + int2string(rightmark);
            else
                glueStrings[kmer].first  = seq + int2string(leftmark) + int2string(rightmark);
        }
        else
        {
            if (left)
                glueStrings[reversecomplement(kmer)].first = reversecomplement(seq) + int2string(rightmark) + int2string(leftmark);
            else
                glueStrings[reversecomplement(kmer)].second = reversecomplement(seq) + int2string(rightmark) + int2string(leftmark);
        }
    }

    void Glue::insert(string seq, bool leftmark, bool rightmark)
    {
        //if (leftmark)
        {
            // don't yet know how to index the hash table with kmers
            //ModelCanon::Kmer kmer = model.codeSeed(seq.substr(0,kmerSize).c_str(), Data::ASCII); 
            string kmer = seq.substr(0,kmerSize);
            aux_addhash(seq, kmer, leftmark, rightmark, true);
        }
        //if (rightmark)
        if (seq.size() > kmerSize)
        {
            //ModelCanon::Kmer kmer = model.codeSeed(seq.substr(seq.size() - kmerSize,kmerSize).c_str(), Data::ASCII); 
            string kmer = seq.substr(seq.size() - kmerSize, kmerSize);
            aux_addhash(seq, kmer, leftmark, rightmark, false);
        }
    }

    string Glue::seq_only(string seq)
    {
        return seq.substr(0,seq.size()-2);
    }

    void Glue::remove(string seq)
    {
        // extremely naive for now
            string leftkmer = seq.substr(0,kmerSize);
            string rightkmer = seq.substr(seq.size() - kmerSize - 2, kmerSize);

            string min = std::min(leftkmer, reversecomplement(leftkmer));
            if (seq_only(glueStrings[min].first) == seq_only(seq) \
                    || seq_only(glueStrings[min].first) == reversecomplement(seq_only(seq)))
            {
                glueStrings[min].first = glueStrings[min].second;
                glueStrings[min].second = "";
            }
            if (seq_only(glueStrings[min].second) == seq_only(seq) \
                    || seq_only(glueStrings[min].second) == reversecomplement(seq_only(seq)))
            {
                glueStrings[min].second = "";
            }

            min = std::min(rightkmer, reversecomplement(rightkmer));
            if (seq_only(glueStrings[min].first) == seq_only(seq) \
                    || seq_only(glueStrings[min].first) == reversecomplement(seq_only(seq)))
            {
                glueStrings[min].first = glueStrings[min].second;
                glueStrings[min].second = "";
            }
            if (seq_only(glueStrings[min].second) == seq_only(seq) \
                    || seq_only(glueStrings[min].second) == reversecomplement(seq_only(seq)))
            {
                glueStrings[min].second = "";
            }

    }

    void Glue::output(string seq)
    {
        Sequence s (Data::ASCII);
        s.getData().setRef ((char*)seq.c_str(), seq.size());
        out.insert(s);
    }

    void Glue::glue()
    {
        size_t k = kmerSize;

        // mostly copypaste of code by Colleen
        for (auto it_gS = glueStrings.begin(); it_gS != glueStrings.end(); it_gS++)
        {
            string query = (it_gS->second).first;   // "query" is always the first node, should be left side of final node
            string match = (it_gS->second).second;

            if (query == "" || match == "")
                continue; // not yet ready to process this

            assert((query.compare("-1") == 0 && match.compare("-1") == 0) || (query.compare("-1") != 0 && match.compare("-1") != 0)); // should never have just one string of pair -1
            // fixme: for some reason assert does not work in this project 
            if (query != "-1" && match=="-1") {printf("assert failed!\n"); exit(1);}

            // already output
            if (query == "-1")
            {
                continue;
            }

            string queryNode = query.substr(0, query.size()-2); // remove leftmark and rightmark
            string rightKmer = queryNode.substr(queryNode.size() - k, k);                       //last k characters of the query node
            string leftKmer = queryNode.substr(0, k);                                           //first k characters of the query node
            string queryLeftMark = query.substr(query.size() - 2, 1);
            string queryRightMark = query.substr(query.size() - 1, 1);

            if (debug)
                cout << "query is: " << queryNode << queryLeftMark << queryRightMark << "\n";
            if (debug)
                cout << "rightKmer is: " << rightKmer << "\n";

            string matchNode = match.substr(0, match.size()-2);
            if (debug)
                cout << "matchnode: '" <<  matchNode << "' match: " << match << std::endl;

            string matchLeftKmer = matchNode.substr(0, k);
            string matchRightKmer = matchNode.substr(matchNode.size() - k, k);
            string matchLeftMark = match.substr(match.size() - 2, 1);
            string matchRightMark = match.substr(match.size() - 1, 1);

            if (matchLeftMark != "1" || rightKmer.compare(matchLeftKmer) != 0) // try swapping
            {
                std::swap(query, match);

                queryNode = query.substr(0, query.size()-2); // remove leftmark and rightmark
                rightKmer = queryNode.substr(queryNode.size() - k, k);                       //last k characters of the query node
                leftKmer = queryNode.substr(0, k);                                           //first k characters of the query node
                queryLeftMark = query.substr(query.size() - 2, 1);
                queryRightMark = query.substr(query.size() - 1, 1);

                matchNode = match.substr(0, match.size()-2);
                matchLeftKmer = matchNode.substr(0, k);
                matchRightKmer = matchNode.substr(matchNode.size() - k, k);
                matchLeftMark = match.substr(match.size() - 2, 1);
                matchRightMark = match.substr(match.size() - 1, 1);

            }

            if (matchLeftMark != "1")
                continue; // nothing to do here

            if (rightKmer.compare(matchLeftKmer) == 0)
            {
                string node = queryNode + matchNode.substr(k, matchNode.size() - k);
                it_gS->second.first = "-1";     // these are done being glued
                it_gS->second.second = "-1";

                remove(query); remove(match); // scrub the table

                if (queryLeftMark == "0" && matchRightMark == "0")
                {
                    if (debug)
                        cout << "Case 1 match and output: " << matchNode << " " << matchLeftMark << matchRightMark << "\n";
                    output(node);
                }
                else    // at least one of left and right needs to be updated
                {
                    insert(node, queryLeftMark == "1", matchRightMark == "1");
                }
            }
            else
            {
                cout << "uh oh, not matching" << endl;
                cout << "  query: " << debug_highlight(query,rightKmer) <<  "\n";
                cout << "  match: " << debug_highlight(match,matchLeftKmer) << std::endl;
            }
        }

        for (auto it = glueStrings.begin(); it != glueStrings.end();)
        {
            if ((it->second).first == "-1" && (it->second).second == "-1")
            {
#ifdef SPARSEHASH
                glueStrings.erase(it);
#else
                it = glueStrings.erase(it);
#endif
            }
            else
                it ++;
        }
#ifdef SPARSEHASH
        glueStrings.resize(0); // effectively remove erased entries from memory
#endif
        cout << "new size: " << glueStrings.size() <<endl;


        /*
        cout << "what remains in the glue\n";
        for (auto it = glueStrings.begin(); it != glueStrings.end(); it++)
        {
            cout << debug_highlight((it->second).first,it->first) << " - " << debug_highlight((it->second).second,it->first) << "\n";
            {
                string seq = it->first;
                Model modelK1(kmerSize-1, minSize);
                Model::Kmer kmmerBegin=modelK1.codeSeed(seq.substr(0,kmerSize-1).c_str(),Data::ASCII);
                size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
                Model::Kmer kmmerEnd=modelK1.codeSeed(seq.substr(seq.size()-kmerSize+1,kmerSize-1).c_str(),Data::ASCII);
                size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));
            }
        }
        */
    }

void put_into_glue(string seq, size_t minimizer, BankFasta &out, Glue &glue, Model& modelK1)
{

    Model::Kmer kmmerBegin=modelK1.codeSeed(seq.substr(0,kmerSize-1).c_str(),Data::ASCII);
    size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
    Model::Kmer kmmerEnd=modelK1.codeSeed(seq.substr(seq.size()-kmerSize+1,kmerSize-1).c_str(),Data::ASCII);
    size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));
    
    size_t i(minimizer/numBucket);
    size_t j(minimizer%numBucket); // following notation from the ismb2015 paper
   
    bool leftmark  = minimizer != leftMin;
    bool rightmark = minimizer != rightMin; 

    //size_t b_pre_major = leftMin/numBucket;
    //size_t b_suf_major = rightMin/numBucket;

    if((!leftmark) && (!rightmark))
    { 
        glue.output(seq);
    } 
    else
    { /* put into glue */
        glue.insert(seq, leftmark, rightmark);
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
    bool parallel = (threads > 1);  

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
            if (parallel)
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

    auto start_buckets=chrono::system_clock::now();

    /**FOREACH SUPERBUCKET **/
    for(uint i(0);i<numBucket;++i){
        cout<<'-'<<flush;
        if(superBuckets[i]!=NULL){
            auto start_flush_t=omp_get_wtime();
            superBuckets[i]->flush();
            auto end_flush_t=omp_get_wtime();
#pragma omp atomic
            global_wtime_flush_sb += omp_diff_wtime_bcalm(start_flush_t, end_flush_t);
            if(superBuckets[i]->getSize()>0){
                auto start_createbucket_t=omp_get_wtime();
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
                    
                    if (parallel)
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
                }
                auto end_createbucket_t=omp_get_wtime();
#pragma omp atomic
                global_wtime_create_buckets += omp_diff_wtime_bcalm(start_createbucket_t, end_createbucket_t);

                auto start_foreach_bucket_t=omp_get_wtime();

                /**FOREACH BUCKET **/
                for(uint j(0);j<numBucket;++j){
                    if(Buckets[j]!=NULL){
                        Buckets[j]->flush();
                        if(Buckets[j]->getSize()>0){
                            //~ graph1 g(kmerSize);
                            
                            /* add nodes to graph */
                            auto start_nodes_t=omp_get_wtime();
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
                                
                                g.addvertex(tmp);
                                g.addleftmin(leftMin);
                                g.addrightmin(rightMin);
                            }
                            auto end_nodes_t=omp_get_wtime();
#pragma omp atomic
                            global_wtime_add_nodes += omp_diff_wtime_bcalm(start_nodes_t, end_nodes_t);

                            /* compact graph*/
                            auto start_dbg=omp_get_wtime();
                            g.debruijn();
                            //~ g.compressh(actualMinimizer);
                            g.compress2();
                            auto end_dbg=omp_get_wtime();
#pragma omp atomic
                            global_wtime_compactions += omp_diff_wtime_bcalm(start_dbg, end_dbg);

                            /* distribute nodes (to other buckets, or output, or glue) */
                            auto start_cdistribution_t=omp_get_wtime(); 
                            for(uint32_t i(1);i<g.unitigs.size();++i){
                                if(g.unitigs[i].size()!=0){
                                    if (parallel)
                                        put_into_glue(g.unitigs[i], actualMinimizer, out, glue, modelK1);
                                    else /* write to another bucket or output */
                                        goodplace(g.unitigs[i],actualMinimizer,out,superBuckets,Buckets,modelK1);

                                }
                            }
                            auto end_cdistribution_t=omp_get_wtime();
#pragma omp atomic
                            global_wtime_cdistribution += omp_diff_wtime_bcalm(start_cdistribution_t, end_cdistribution_t);

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
                            delete(Buckets[j]);
                            remove(("B"+to_string(j)).c_str());
                        }
                    }
                }
                delete(superBuckets[i]);
                remove(("SB"+to_string(i)).c_str());

                auto end_foreach_bucket_t=omp_get_wtime();
        #pragma omp atomic
                global_wtime_foreach_bucket += omp_diff_wtime_bcalm(start_foreach_bucket_t, end_foreach_bucket_t);

                // do the gluing at the end of each superbucket
                if (parallel) {
                    auto start_glue_t=omp_get_wtime();

					glue.updateMemStats();
					glue.glue();

                    auto end_glue_t=omp_get_wtime();
#pragma omp atomic
                    global_wtime_glue += omp_diff_wtime_bcalm(start_glue_t, end_glue_t);
 
				}
            }
        }
    }
	
	glue.printMemStats();

    /* printing some timing stats */
    auto end_t=chrono::system_clock::now();
    cout<<"Buckets compaction and gluing: "<<chrono::duration_cast<chrono::nanoseconds>(end_t - start_buckets).count() / unit<<" secs"<<endl;
    cout<<"Within that, \n";
    cout <<"     creating buckets from superbuckets: "<< global_wtime_create_buckets / unit <<" secs"<<endl;
    cout <<"     bucket compaction (except gluing): "<< global_wtime_foreach_bucket / unit <<" secs" <<endl;
    cout <<"     within that, \n";
    cout << "                 flushing superbuckets: "<< global_wtime_flush_sb / unit <<" secs" <<endl;
    cout << "                 adding nodes to subgraphs: "<< global_wtime_add_nodes / unit <<" secs" <<endl;
    cout <<"                  subgraphs constructions and compactions: "<< global_wtime_compactions / unit <<" secs"<<endl;
    cout <<"                  compacted nodes redistribution: "<< global_wtime_cdistribution / unit <<" secs"<<endl;
    cout <<"     glueing "<< global_wtime_glue / unit <<" secs"<<endl;
    double sum = global_wtime_glue + global_wtime_cdistribution + global_wtime_compactions + global_wtime_add_nodes + global_wtime_flush_sb + global_wtime_create_buckets;
    cout<<"Sum of the above fine-grained timings: "<< sum / unit <<" secs"<<endl;
    cout<<"BCALM total wallclock (incl kmer counting): "<<chrono::duration_cast<chrono::nanoseconds>(end_t-start_t).count() / unit <<" secs"<<endl;
    cout<<"Discrepancy between sum of fine-grained timings and total wallclock of buckets compactions step: "<< (chrono::duration_cast<chrono::nanoseconds>(end_t-start_buckets).count() - sum ) / unit <<" secs"<<endl;
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
