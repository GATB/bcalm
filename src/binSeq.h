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
// #define kmer uint64_t
//if 32<k<=64
#define kmer __uint128_t


static uint kmerBinseq;


using namespace std;


class binSeq{
public:
	vector<uint8_t> vect;
	bool isInt;

	string str();
	void str2(string& res);
	binSeq sub();
	binSeq sub(uint begin,uint size);
	binSeq getBegin(uint size);
	binSeq getEnd(uint size);
	kmer getEndInt();
	kmer getBeginInt();
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


void initBinSeq(uint k);


void testBinSeq();


#endif /* defined(__PBMOG__binSeq__) */
