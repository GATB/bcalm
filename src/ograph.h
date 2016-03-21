#ifndef OGRAPH
#define OGRAPH

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cstdlib>
#include "binSeq.h"


using namespace std;


struct comparator{bool operator()(const kmerIndice& a , const kmerIndice& b) { return a.kmmer < b.kmmer; }};


class graph3{
	public:
		uint k,indiceUnitigs,nbElement,minimizer,minsize;
		string* unitigs;
		vector<kmerIndice> left;
		vector<kmerIndice> right;
		void addvertex(string& str);
        void addtuple(tuple<string,uint,uint>& tuple);
		void addleftmin(unsigned int mini);
		void addrightmin(unsigned int mini);
		void debruijn();
        void debruijn2();
        void compaction2(uint iL, uint iR);
		void compress();
        kmer end2int128rc(const string& str);
        kmer end2int128(const string& str);
        kmer beg2int128rc(const string& str);
        kmer beg2int128(const string& str);
        kmer rcb(kmer min);
		void compaction(uint iR, uint iL);
		uint size();
        bool output(uint i);
        bool clear();

		graph3(uint ka, uint min,uint size, uint nb){
            indiceUnitigs=0;
			minsize=size;
			k=ka;
			minimizer=min;
            nbElement=nb;
            unitigs=new string [nbElement];
            left.reserve(nbElement);
    	    right.reserve(nbElement);
		}
};

uint chartoint(char c);
string reverseinplace(string& str);
void reverseinplace2(string& str);
bool isNumber(char c);

#endif
