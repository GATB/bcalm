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


class graph3{
	public:
		uint k,indiceUnitigs,nbElement,minimizer,minsize;
        // uint ;
		// uint32_t minimizer;
		// uint32_t minsize;
        // uint32_t nbElement;
        // uint lol;
        // bool found;
        // uint u1,u2,u3,u4;
        // kmerIndice ki;
        // kmer kmer1,kmer2;
        // ,beg2,beg1,begrc2,end1,end2,endrc2, resrcb,offsetrcb, resBeg,resEnd;
		string* unitigs;
        // string unitigL,unitigR;
		vector<kmerIndice> left;
		vector<kmerIndice> right;
		// vector<kmer2Indice> compactions;
		void addvertex(string& str);
        void addtuple(tuple<string,uint,uint>& tuple);
		void addleftmin(unsigned int mini);
		void addrightmin(unsigned int mini);
		void debruijn();
        void debruijn2();
        void compaction2(uint iL, uint iR);
        // void compaction(const uint32_t& iL, const uint32_t& iR);
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
            // u1=u2=u3=u4=indiceUnitigs=indiceRight=indiceLeft=indiceAdd=0;
            indiceUnitigs=0;
			minsize=size;
			k=ka;
            // found=false;
			minimizer=min;
            nbElement=nb;
            unitigs=new string [nbElement];
            left.reserve(nbElement);
    	    right.reserve(nbElement);
            // compactions.reserve(nbElement);
		}
};

uint chartoint(char c);
string reverseinplace(string& str);
void reverseinplace2(string& str);
bool isNumber(char c);

#endif
