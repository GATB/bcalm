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
class GlueEntry {
	public:

		GlueEntry() : seq(""), lmark(false), rmark(false), kmerSize(-1) { };
		GlueEntry(string raw, size_t _kmerSize) : kmerSize(_kmerSize) {
			seq = raw.substr(0, raw.length() - 2);
			lmark = raw.at(raw.length() - 2);
			rmark = raw.at(raw.length() - 1);
		}
		GlueEntry(string _seq, bool _lmark, bool _rmark, size_t _kmerSize) : seq(_seq), lmark(_lmark), rmark(_rmark), kmerSize(_kmerSize) { };

		string getSeq()   { return seq;   };
		bool getLmark()   { return lmark; }; 
		bool getRmark()   { return rmark; };
		string getLkmer() { return seq.substr(0, kmerSize);  };
		string getRkmer() { return seq.substr(seq.length() - kmerSize, kmerSize);  };
		bool isEmpty()    { return seq == ""; };

		string getRaw() {
			string retval;
			retval = seq;
			retval.push_back(lmark);
			retval.push_back(rmark);
			return retval;
		}

		void revComp() {
			seq = reversecomplement(seq);
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

		void output(string seq);
		void insert_aux(GlueEntry newEntry, string key, GlueEntry & glueResult); 
		void glueSingleEntry(GlueEntry query, GlueEntry match, string key, GlueEntry & glueResult);
		bool check_if_empty(GlueEntry newEntry, string key);

		//timing code
		std::chrono::system_clock::time_point startTime;
		std::chrono::seconds totalTime;
		int timerReferenceCount = 0;
		void startTimer() { 
			if (timerReferenceCount++ == 0) startTime =  chrono::system_clock::now(); 
		};
		void stopTimer();
};

