//
//  binSeq.h
//  PBMOG
//
//  Created by malfoy on 27/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#ifndef __PBMOG__binSeq__
#define __PBMOG__binSeq__

#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <unordered_map>

typedef unsigned int  uint;
//if k<=32
// #define kmer __uint128_t
//if 32<k<=64
#define kmer __uint128_t


struct kmerIndice{
	uint32_t indice;
	kmer kmmer;
};


using namespace std;


class binSeq{
public:
	//~ static const unsigned char rc[];
	//~ static const unordered_map<string,unsigned char> string2Byte;
	//~ bool isNumber;
	vector<uint8_t> vect;
	bool isInt;

	string str();
	void str2(string& res);
	binSeq sub(uint begin);
	binSeq sub(uint begin,uint size);
	binSeq getBegin(uint size);
	binSeq getEnd(uint size);
	kmer getEndInt(uint size);
	kmer getBeginInt(uint size);
	kmer getBeginRcInt(uint size);
	kmer getEndRcInt(uint size);
	void reverse();
	binSeq getReverse();
	void add(const binSeq& b);
	void resize();
	void clear();
	uint size();
	uint getNumber();
	kmer getInt();

	explicit binSeq(const string& str);
	binSeq();
	binSeq(const binSeq& bs);
	explicit binSeq(uint);
};


void initBinSeq();


void testBinSeq();


#endif /* defined(__PBMOG__binSeq__) */
