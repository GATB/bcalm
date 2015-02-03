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
	vector<uint8_t> vect;
	bool isNumber;

	string str();
	binSeq sub(size_t begin);
	binSeq sub(size_t begin,size_t size);
	binSeq getBegin(size_t size);
	binSeq getEnd(size_t size);
//	uint64_t getBegin(size_t size);
//	uint64_t getEnd(size_t size);
//	uint64_t getBeginRc(size_t size);
//	uint64_t getEndRc(size_t size);
	void reverse();
	binSeq getReverse();
	void add(const binSeq& b);
	void resize();
	void clear();
	size_t size();
	uint32_t getNumber();
	uint64_t getInt();

	binSeq& operator+=(const binSeq& rhs){
		this->add(rhs);
		return *this;
	}






	binSeq(const string& str);
	binSeq();
	binSeq(const binSeq& bs);
	binSeq(uint32_t);


};


inline bool operator==(const binSeq& lhs, const binSeq& rhs){return lhs.vect==rhs.vect;}
inline bool operator!=(const binSeq& lhs, const binSeq& rhs){return !operator==(lhs,rhs);}
inline binSeq operator+(binSeq lhs, const binSeq& rhs){
	lhs += rhs;
	return lhs;
}



unsigned char char2int(unsigned char c);
unsigned char int2char(unsigned char c);

void testBinSeq();


#endif /* defined(__PBMOG__binSeq__) */
