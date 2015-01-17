#include <assert.h>
#include <gatb/gatb_core.hpp>
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

// glue has threaded components now
#include <thread>
#include <atomic>
#include <mutex>
#include <../../../thirdparty/concurrentqueue.h>

#define TWOBITGLUEHASH 
//#define SPARSEHASH 
#ifdef SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;
#endif 

using namespace std;

string revcomp (string &s);
//string reversecomplement(const string& dna); // code taken from TakeABreak
string debug_highlight(string s, string motif);

#ifdef TWOBITGLUEHASH
typedef vector<uint8_t> RawEntry;
#else
typedef string RawEntry;
#endif

//class to represent an entry in GlueStorage
class GlueEntry {
	public:

		GlueEntry() : seq(""), lmark(false), rmark(false), kmerSize(-1) { };
		GlueEntry(string _seq, bool _lmark, bool _rmark, size_t _kmerSize) : seq(_seq), lmark(_lmark), rmark(_rmark), kmerSize(_kmerSize) { };
		GlueEntry(RawEntry raw, size_t _kmerSize);

		string getSeq()   { return seq;   };
		bool getLmark()   { return lmark; }; 
		bool getRmark()   { return rmark; };
		string getLkmer() { return seq.substr(0, kmerSize);  };
		string getRkmer() { return seq.substr(seq.length() - kmerSize, kmerSize);  };
		bool isEmpty()    { return seq == ""; };

		RawEntry getRaw();

		void revComp() {
			seq = revcomp(seq);
			std::swap(lmark, rmark);
		}

	private:
		string seq;
		bool lmark;
		bool rmark;
		size_t kmerSize;
};

string tostring(const GlueEntry e, string key);

/*
class TagSystem {
	public:
		typedef long int TagIndex;

		TagIndex createTag(GlueEntry e);

	private:

};
*/




class GlueStorage {
	public:
		int kmerSize;

		GlueStorage(int _kmerSize) : kmerSize(_kmerSize) { 
#ifdef SPARSEHASH
			glueMap.set_deleted_key("");
#endif
		}
		bool find (string key, GlueEntry & e);
		//bool getFirst(GlueEntry & e1, GlueEntry & e2);
		//bool getNext(GlueEntry & e1, GlueEntry & e2); 
		//void insertAtCurIt(GlueEntry e1, GlueEntry e2) { insertAtIt(curIt, e1, e2); }
		void insertAtKey(string key, GlueEntry e);
		void insertAfterFind(GlueEntry e) { insertAtIt(findIt, e); } 
		void cleanup();
		void updateMemStats();
		void printMemStats();
        unsigned long glueMapSize();
        unsigned long glueMapSizeNonEmpty();

		string dump();
		string dump(string key, bool dumpkey = true);
        
        unsigned long nbNonEmpty;

	private:
#ifdef SPARSEHASH
		typedef  sparse_hash_map<string, RawEntry> GlueMap; 
#else
		typedef unordered_map<string, RawEntry> GlueMap;
#endif
		GlueMap glueMap;
		GlueMap::iterator curIt; //iterator state used by getFirst and getNext
		GlueMap::iterator findIt; //iterator state used by find and insertAfterFind

		//used for tracking memory usage
		size_t maxEntries = 0;
		size_t totEntries = 0;
		int numDataPoints = 0;
		size_t maxSize = 0;
		size_t totSize = 0;

		bool derefIt (GlueMap::const_iterator it, GlueEntry & e);
		void insertAtIt(GlueMap::iterator it, GlueEntry e);
};


class Glue;


class GlueCommander
{
    public:
        // to compute minimizers
        const static size_t SPAN = KSIZE_2; // TODO: should be shared with bcalm_1.cpp
        typedef Kmer<SPAN>::ModelCanonical ModelCanon;
        typedef Kmer<SPAN>::ModelMinimizer <ModelCanon> Model;


    void insert(GlueEntry &e);
    bool insert_aux(GlueEntry newEntry, string key);
    void cleanup_force();
    void cleanup_threaded();
    void updateMemStats();
    void dump();
    void printMemStats();
    void stop();
	void output(string seq);
    void spawn_threads();
    int which_queue(size_t minimizer);
    unsigned long queues_size(bool silent=false);

	GlueCommander(size_t _kmerSize, BankFasta *out, int nb_glues, Model *model);


    private:
        std::mutex commander_mutex;
        BankFasta *out;
        Model *model;
        int nb_glues;
        std::vector<Glue*> glues;
        std::vector<moodycamel::ConcurrentQueue<pair<GlueEntry, string> >>  insert_aux_queues;
        std::vector<bool> keep_glueing; 
        std::vector<std::thread*> glue_threads; 
};

class Glue
{
    public:
		GlueStorage glueStorage; //this should really be treated as private. It is only public to allow calling updateMemStats and such

		Glue(size_t _kmerSize, GlueCommander *commander) : kmerSize(_kmerSize), glueStorage(_kmerSize), commander(commander), nbGlueInserts(0), needs_cleanup(false) {
#ifdef TWOBITGLUEHASH 
			cout << "Glue: using TWOBITBLUEHASH.\n";
#else
			cout << "Glue: using a string hash.\n";
#endif
#ifdef SPARSEHASH
			cout << "Glue: using SPARSEHASH.\n";
#else
			cout << "Glue: using unordered_map.\n";
#endif
		}
		void insert(GlueEntry newEntry, bool process = false);
		void glue();
		
		unsigned long getTime() { 
			if (timerReferenceCount != 0) {
				cout << "error in Glue::getTime(): timerReference count is " << timerReferenceCount << endl;
				exit(1);
			}
			return totalTime.count();
		};

		bool insert_aux(GlueEntry newEntry, string key, GlueEntry & glueResult, bool onlyCheckIfEmpty = false);
        
        std::atomic<bool> needs_cleanup;

	private:
		int kmerSize;
        GlueCommander *commander;

		void glueSingleEntry(GlueEntry query, GlueEntry match, string key, GlueEntry & glueResult);
		bool check_if_empty(GlueEntry newEntry, string key);

		//timing code
		std::chrono::system_clock::time_point startTime;
		std::chrono::microseconds totalTime;
		int timerReferenceCount = 0;
		void startTimer() { 
			if (timerReferenceCount++ == 0) startTime =  chrono::system_clock::now(); 
		};
		void stopTimer();

        unsigned long nbGlueInserts;
};

