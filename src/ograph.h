#ifndef OGRAPH
#define OGRAPH

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cstdlib>
#include "binSeq.h"

#ifdef SPARSE_HASH
    #include <sparsehash/sparse_hash_map>
    #include <sparsehash/dense_hash_map>
    #include <sparsehash/dense_hash_set>

    #define unordered_map sparse_hash_map
    #define unordered_map dense_hash_map
    #define unordered_set dense_hash_set
#else
    #include <unordered_map>
    #include <unordered_set>

#endif /* SPARSE_HASH */

using namespace std;

void create_hash_function_from_m_mers(int m);
void count_m_mers(string str, int m, int k);
void init_m_mers_table(int m);

typedef unordered_map<string,int> HashMap;


HashMap build_hash_map(int len);

int shash(const string& s, int& previous_shash, unsigned int start_pos = 0, int length = -1);

string inverse_shash (int num, int len);

int minimiserrc(const string &node,const int &minimisersize);

int minbutbiggerthan(int leftmin, int rightmin, const string &namebucket);

string reversecompletment(const string& str );

bool adjacent (const string& node1,const  string& node2,int k);

string readn(ifstream *file,uint64_t n);

string minimalsub(const string &w, const int &p,const int &k);

string minimalsub2(const string &w, const int &p,const int &k);

class neighbour
{
	public:
		array<pair<uint64_t,unsigned char>,8> list;
		//~ neighbour()
		//~ {
			//~ for(int )
			//~ list[i]=make_pair(0,0);
		//~ }
		uint64_t nbtype(unsigned char c);
		uint64_t gtype(unsigned char c);

		void add(uint64_t p,unsigned char b);
		unsigned char remove(uint64_t v);
		unsigned char removep(uint64_t v,unsigned char c);
		unsigned char removetype(unsigned char c);

};


class graph2
{
	public:
		int k;
		int minimizer;
		int minsize;
		uint64_t weight;
		//~ ofstream outFile;
		vector<string> unitigs;
		vector<bool> leftmins;
		vector<bool> rightmins;
		unordered_map<uint64_t,uint32_t> left2unitig;
		unordered_map<uint64_t,uint32_t> right2unitig;
		//~ dense_hash_set<uint64_t> multiple;
		//~ unordered_set<uint64_t> multiple;

		void addvertex(const string& str);
		void addleftmin(int mini);
		void addrightmin(int mini);
		void debruijn();
		void compress();
		void print();
		void chainCompaction(uint32_t i,string unitig, uint32_t next);
		void chainCompaction2(uint32_t i,string unitig, uint32_t next);
		uint32_t leftUnique(uint64_t);
		uint32_t rightUnique(uint64_t);
		uint32_t goBeg(uint32_t i);
		uint32_t goEnd(uint32_t i);
		uint32_t size();

		graph2(int ka, int min,int size)
		{
			minsize=size;
			k=ka;
			minimizer=min;
			//~ multiple.set_empty_key(-1);
//			left2unitig.set_empty_key(-1);
//			right2unitig.set_empty_key(-1);
			unitigs.push_back("");
			leftmins.push_back(-1);
			rightmins.push_back(-1);
		}
};


struct kmerIndice{
	uint32_t indice;
	__uint128_t kmmer;

};

struct kmer2Indice{
	uint32_t indiceL;
	uint32_t indiceR;
};



class graph3
{
	public:
		uint k;
        uint indiceAdd,indiceLeft,indiceRight,indiceUnitigs;
		uint32_t minimizer;
		uint32_t minsize;
        uint32_t nbElement;
        kmerIndice ki;
        // __uint128_t kmer1,kmer2;
        // ,beg2,beg1,begrc2,end1,end2,endrc2, resrcb,offsetrcb, resBeg,resEnd;
		string* unitigs;
		vector<kmerIndice> left;
		vector<kmerIndice> right;
		// vector<kmer2Indice> compactions;
		void addvertex(string& str);
        void addtuple(tuple<string,uint32_t,uint32_t>& tuple);
		void addleftmin(unsigned int mini);
		void addrightmin(unsigned int mini);
		void debruijn();
        void debruijn2();
        void compaction2(uint32_t iL, uint32_t iR);
        // void compaction(const uint32_t& iL, const uint32_t& iR);
		void compress();
        __uint128_t end2int128rc(const string& str);
        __uint128_t end2int128(const string& str);
        __uint128_t beg2int128rc(const string& str);
        __uint128_t beg2int128(const string& str);
        __uint128_t rcb(__uint128_t min);
		void compaction(uint32_t iR, uint32_t iL);
		uint32_t size();
        bool output(uint i);
        bool clear();

		graph3(uint32_t ka, uint32_t min,uint32_t size, uint nb){
            indiceUnitigs=indiceRight=indiceLeft=indiceAdd=0;
			minsize=size;
			k=ka;
			minimizer=min;
            nbElement=nb;
            unitigs =new string [nbElement];
            left.reserve(nbElement);
    	    right.reserve(nbElement);
            // compactions.reserve(nbElement);
		}
};

class graph4
{
	public:
		//~ uint32_t a,b,c,d;
		uint32_t k;
		uint32_t minimizer;
		uint32_t minsize;
		__uint128_t maximum;
		//~ uint32_t nbKmer;
		vector<binSeq> unitigs;
		vector<bool> leftmins;
		vector<bool> rightmins;
		vector<bool> isNumber;
		vector<kmerIndice> left;
		vector<kmerIndice> right;
		//~ vector<kmer2Indice> compactions;

		void addvertex(const string& str);
		void addleftmin(unsigned int mini);
		void addrightmin(unsigned int mini);
		void debruijn();

		void compress();
		void compaction(uint32_t iR, uint32_t iL);
		uint32_t size();

		graph4(uint32_t ka, uint32_t min,uint32_t size){
			//~ a=b=c=d=0;
			//~ nbKmer=kmerInGraph;
			maximum=0;
			minsize=size;
			k=ka;
			minimizer=min;
			//~ unitigs.reserve(nbKmer);
			//~ leftmins.reserve(nbKmer);
			//~ rightmins.reserve(nbKmer);
			//~ left.reserve(nbKmer);
			//~ right.reserve(nbKmer);
			//~ isNumber.assign(nbKmer,false);
		}
};

//~ class graph5
//~ {
	//~ public:
		//~ uint32_t k;
		//~ uint32_t minimizer;
		//~ uint32_t minsize;
		//~ uint64_t maxR,maxL;
		//~ vector<binSeq2> unitigs;
		//~ vector<bool> leftmins;
		//~ vector<bool> rightmins;
		//~ vector<kmer2Indice> compactions;
//~
		//~ void addvertex(const string& str);
		//~ void addleftmin(int mini);
		//~ void addrightmin(int mini);
		//~ void debruijn();
		//~ void compress();
		//~ void compaction(uint32_t iR, uint32_t iL);
		//~ uint32_t size();
//~
		//~ graph4(uint32_t ka, uint32_t min,uint32_t size){
			//~ maxR=maxL=0;
			//~ minsize=size;
			//~ k=ka;
			//~ minimizer=min;
		//~ }
//~ };

#endif
