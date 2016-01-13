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
#include <unordered_map>

using namespace std;


class binSeq {
public:
	//~ static const unsigned char rc[];
	//~ static const unordered_map<string,unsigned char> string2Byte;
	//~ bool isNumber;
	vector<uint8_t> vect;

	bool isInt;
	string str();
	void str2(string& res);
	binSeq sub(uint8_t begin);
	binSeq sub(size_t begin,size_t size);
	binSeq getBegin(uint8_t size);
	binSeq getEnd(uint8_t size);
	__uint128_t getEndInt(uint8_t size);
	__uint128_t getBeginInt(uint8_t size);
	__uint128_t getBeginRcInt(uint8_t size);
	__uint128_t getEndRcInt(uint8_t size);
	void reverse();
	binSeq getReverse();
	void add(const binSeq& b);
	void resize();
	void clear();
	size_t size();
	uint32_t getNumber();
	__uint128_t getInt();

	explicit binSeq(const string& str);
	binSeq();
	binSeq(const binSeq& bs);
	explicit binSeq(uint32_t);
};


//~ inline bool operator==(const binSeq& lhs, const binSeq& rhs){return lhs.vect==rhs.vect;}
//~ inline bool operator!=(const binSeq& lhs, const binSeq& rhs){return !operator==(lhs,rhs);}
//~ inline binSeq operator+(binSeq lhs, const binSeq& rhs){
	//~ lhs += rhs;
	//~ return lhs;
//~ }


void initBinSeq();


void testBinSeq();


#endif /* defined(__PBMOG__binSeq__) */
