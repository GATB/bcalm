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

//#define SPARSEHASH 
#ifdef SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;
#endif 

using namespace std;

string reversecomplement(const string& dna); // code taken from TakeABreak
string debug_highlight(string s, string motif);

//class to represent an entry in GlueStorage
//it is meant to be used only when an entry is extracted from storage. 
//For representation in storage, we use GlueEntryCompact*
struct GlueEntry {
	string seq = "";
	bool lmark = 0;
	bool rmark = 0;
	string lkmer = "";
	string rkmer = "";
};

string tostring(const GlueEntry e, string key);

// class to support compact representation of a GlueEntry, as should be stored in memory
// For now, this is very naive, and just stores it a string
class GlueEntryCompactNaive{
	public:
		GlueEntryCompactNaive(string _raw, size_t _kmerSize) : raw(_raw), kmerSize(_kmerSize) { }
		GlueEntryCompactNaive(GlueEntry e, size_t _kmerSize) : kmerSize(_kmerSize) {
			raw = e.seq;
			raw.push_back(e.lmark); 
			raw.push_back(e.rmark);
		}
		string getRaw() { return raw; }
		GlueEntry getEntry();

	private:
		string raw;
		size_t kmerSize;
};



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

		string dump();
		string dump(string key, bool dumpkey = true);

	private:
#ifdef SPARSEHASH
		typedef  sparse_hash_map<string, string> GlueMap; 
#else
		typedef unordered_map<string, string> GlueMap;
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


class Glue
{
    public:
		GlueStorage glueStorage; //this should really be treated as private. It is only public to allow calling updateMemStats and such

		Glue(size_t _kmerSize, BankFasta &out) : kmerSize(_kmerSize), out(out), glueStorage(_kmerSize) { }
		void insert(GlueEntry newEntry, bool process = false);
		void glue();
		
		unsigned long getTime() { 
			if (timerReferenceCount != 0) {
				cout << "error in Glue::getTime(): timerReference count is " << timerReferenceCount << endl;
				exit(1);
			}
			return totalTime.count();
		};

	private:
		BankFasta out;
		int kmerSize;

		std::chrono::system_clock::time_point startTime;
		std::chrono::seconds totalTime;
		int timerReferenceCount = 0;

		void output(string seq);
		void insert_aux(GlueEntry newEntry, string key, GlueEntry & glueResult); 
		bool glueSingleEntry(GlueEntry query, GlueEntry match, string key, GlueEntry & glueResult);
		void startTimer() { 
			if (timerReferenceCount++ == 0) startTime =  chrono::system_clock::now(); 
		};
		void stopTimer();
};

