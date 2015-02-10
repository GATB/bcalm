//
//  binSeq.cpp
//  PBMOG
//
//  Created by malfoy on 27/01/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "binSeq.h"
#include <string>
#include <iostream>
#include <algorithm>




static const uint8_t rc[]={

	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,

	0b01000011,0b01000010,0b01000001,0b01000000,0b01000000,0b01000000,0b01000000,0b01000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,

	0b10001111,0b10001011,0b10000111,0b10000011,0b10001110,0b10001010,0b10000110,0b10000010,
	0b10001101,0b10001001,0b10000101,0b10000001,0b10001100,0b10001000,0b10000100,0b10000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
	0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,

	0b11111111,0b11101111,0b11011111,0b11001111,0b11111011,0b11101011,0b11011011,0b11001011,
	0b11110111,0b11100111,0b11010111,0b11000111,0b11110011,0b11100011,0b11010011,0b11000011,
	0b11111110,0b11101110,0b11011110,0b11001110,0b11111010,0b11101010,0b11011010,0b11001010,
	0b11110110,0b11100110,0b11010110,0b11000110,0b11110010,0b11100010,0b11010010,0b11000010,
	0b11111101,0b11101101,0b11011101,0b11001101,0b11111001,0b11101001,0b11011001,0b11001001,
	0b11110101,0b11100101,0b11010101,0b11000101,0b11110001,0b11100001,0b11010001,0b11000001,
	0b11111100,0b11101100,0b11011100,0b11001100,0b11111000,0b11101000,0b11011000,0b11001000,
	0b11110100,0b11100100,0b11010100,0b11000100,0b11110000,0b11100000,0b11010000,0b11000000,

};


static uint8_t char2int1[100];

static uint8_t char2int[255];

static uint8_t char2int2[255];

uint8_t functionwihtoutcoll(uint8_t x,uint8_t y,uint8_t z){
	return( (x&0b00000110)|((y&0b00000110)<<2)|((z&0b00000110)<<4) );
}

uint8_t functionwihtoutcoll(uint8_t x,uint8_t y){
	return(x|(y<<4));
}


void initBinSeq(){
	char2int1[65]=0b01000000;
	char2int1[67]=0b01000001;
	char2int1[71]=0b01000010;
	char2int1[84]=0b01000011;

	char2int2[functionwihtoutcoll(65,65)]=0b10000000;
	char2int2[functionwihtoutcoll(65,67)]=0b10000001;
	char2int2[functionwihtoutcoll(65,71)]=0b10000010;
	char2int2[functionwihtoutcoll(65,84)]=0b10000011;
	char2int2[functionwihtoutcoll(67,65)]=0b10000100;
	char2int2[functionwihtoutcoll(67,67)]=0b10000101;
	char2int2[functionwihtoutcoll(67,71)]=0b10000110;
	char2int2[functionwihtoutcoll(67,84)]=0b10000111;
	char2int2[functionwihtoutcoll(71,65)]=0b10001000;
	char2int2[functionwihtoutcoll(71,67)]=0b10001001;
	char2int2[functionwihtoutcoll(71,71)]=0b10001010;
	char2int2[functionwihtoutcoll(71,84)]=0b10001011;
	char2int2[functionwihtoutcoll(84,65)]=0b10001100;
	char2int2[functionwihtoutcoll(84,67)]=0b10001101;
	char2int2[functionwihtoutcoll(84,71)]=0b10001110;
	char2int2[functionwihtoutcoll(84,84)]=0b10001111;

	char2int[functionwihtoutcoll(65,65,65)]=0b11000000;
	char2int[functionwihtoutcoll(65,65,67)]=0b11000001;
	char2int[functionwihtoutcoll(65,65,71)]=0b11000010;
	char2int[functionwihtoutcoll(65,65,84)]=0b11000011;
	char2int[functionwihtoutcoll(65,67,65)]=0b11000100;
	char2int[functionwihtoutcoll(65,67,67)]=0b11000101;
	char2int[functionwihtoutcoll(65,67,71)]=0b11000110;
	char2int[functionwihtoutcoll(65,67,84)]=0b11000111;
	char2int[functionwihtoutcoll(65,71,65)]=0b11001000;
	char2int[functionwihtoutcoll(65,71,67)]=0b11001001;
	char2int[functionwihtoutcoll(65,71,71)]=0b11001010;
	char2int[functionwihtoutcoll(65,71,84)]=0b11001011;
	char2int[functionwihtoutcoll(65,84,65)]=0b11001100;
	char2int[functionwihtoutcoll(65,84,67)]=0b11001101;
	char2int[functionwihtoutcoll(65,84,71)]=0b11001110;
	char2int[functionwihtoutcoll(65,84,84)]=0b11001111;
	char2int[functionwihtoutcoll(67,65,65)]=0b11010000;
	char2int[functionwihtoutcoll(67,65,67)]=0b11010001;
	char2int[functionwihtoutcoll(67,65,71)]=0b11010010;
	char2int[functionwihtoutcoll(67,65,84)]=0b11010011;
	char2int[functionwihtoutcoll(67,67,65)]=0b11010100;
	char2int[functionwihtoutcoll(67,67,67)]=0b11010101;
	char2int[functionwihtoutcoll(67,67,71)]=0b11010110;
	char2int[functionwihtoutcoll(67,67,84)]=0b11010111;
	char2int[functionwihtoutcoll(67,71,65)]=0b11011000;
	char2int[functionwihtoutcoll(67,71,67)]=0b11011001;
	char2int[functionwihtoutcoll(67,71,71)]=0b11011010;
	char2int[functionwihtoutcoll(67,71,84)]=0b11011011;
	char2int[functionwihtoutcoll(67,84,65)]=0b11011100;
	char2int[functionwihtoutcoll(67,84,67)]=0b11011101;
	char2int[functionwihtoutcoll(67,84,71)]=0b11011110;
	char2int[functionwihtoutcoll(67,84,84)]=0b11011111;
	char2int[functionwihtoutcoll(71,65,65)]=0b11100000;
	char2int[functionwihtoutcoll(71,65,67)]=0b11100001;
	char2int[functionwihtoutcoll(71,65,71)]=0b11100010;
	char2int[functionwihtoutcoll(71,65,84)]=0b11100011;
	char2int[functionwihtoutcoll(71,67,65)]=0b11100100;
	char2int[functionwihtoutcoll(71,67,67)]=0b11100101;
	char2int[functionwihtoutcoll(71,67,71)]=0b11100110;
	char2int[functionwihtoutcoll(71,67,84)]=0b11100111;
	char2int[functionwihtoutcoll(71,71,65)]=0b11101000;
	char2int[functionwihtoutcoll(71,71,67)]=0b11101001;
	char2int[functionwihtoutcoll(71,71,71)]=0b11101010;
	char2int[functionwihtoutcoll(71,71,84)]=0b11101011;
	char2int[functionwihtoutcoll(71,84,65)]=0b11101100;
	char2int[functionwihtoutcoll(71,84,67)]=0b11101101;
	char2int[functionwihtoutcoll(71,84,71)]=0b11101110;
	char2int[functionwihtoutcoll(71,84,84)]=0b11101111;
	char2int[functionwihtoutcoll(84,65,65)]=0b11110000;
	char2int[functionwihtoutcoll(84,65,67)]=0b11110001;
	char2int[functionwihtoutcoll(84,65,71)]=0b11110010;
	char2int[functionwihtoutcoll(84,65,84)]=0b11110011;
	char2int[functionwihtoutcoll(84,67,65)]=0b11110100;
	char2int[functionwihtoutcoll(84,67,67)]=0b11110101;
	char2int[functionwihtoutcoll(84,67,71)]=0b11110110;
	char2int[functionwihtoutcoll(84,67,84)]=0b11110111;
	char2int[functionwihtoutcoll(84,71,65)]=0b11111000;
	char2int[functionwihtoutcoll(84,71,67)]=0b11111001;
	char2int[functionwihtoutcoll(84,71,71)]=0b11111010;
	char2int[functionwihtoutcoll(84,71,84)]=0b11111011;
	char2int[functionwihtoutcoll(84,84,65)]=0b11111100;
	char2int[functionwihtoutcoll(84,84,67)]=0b11111101;
	char2int[functionwihtoutcoll(84,84,71)]=0b11111110;
	char2int[functionwihtoutcoll(84,84,84)]=0b11111111;
}



const unordered_map<string,unsigned >  string2Byte={
	{"A",0b01000000},{"C",0b01000001},{"G",0b01000010},{"T",0b01000011},

	{"AA",0b10000000},{"AC",0b10000001},{"AG",0b10000010},{"AT",0b10000011},
	{"CA",0b10000100},{"CC",0b10000101},{"CG",0b10000110},{"CT",0b10000111},
	{"GA",0b10001000},{"GC",0b10001001},{"GG",0b10001010},{"GT",0b10001011},
	{"TA",0b10001100},{"TC",0b10001101},{"TG",0b10001110},{"TT",0b10001111},

	{"AAA",0b11000000},{"AAC",0b11000001},{"AAG",0b11000010},{"AAT",0b11000011},
	{"ACA",0b11000100},{"ACC",0b11000101},{"ACG",0b11000110},{"ACT",0b11000111},
	{"AGA",0b11001000},{"AGC",0b11001001},{"AGG",0b11001010},{"AGT",0b11001011},
	{"ATA",0b11001100},{"ATC",0b11001101},{"ATG",0b11001110},{"ATT",0b11001111},

	{"CAA",0b11010000},{"CAC",0b11010001},{"CAG",0b11010010},{"CAT",0b11010011},
	{"CCA",0b11010100},{"CCC",0b11010101},{"CCG",0b11010110},{"CCT",0b11010111},
	{"CGA",0b11011000},{"CGC",0b11011001},{"CGG",0b11011010},{"CGT",0b11011011},
	{"CTA",0b11011100},{"CTC",0b11011101},{"CTG",0b11011110},{"CTT",0b11011111},

	{"GAA",0b11100000},{"GAC",0b11100001},{"GAG",0b11100010},{"GAT",0b11100011},
	{"GCA",0b11100100},{"GCC",0b11100101},{"GCG",0b11100110},{"GCT",0b11100111},
	{"GGA",0b11101000},{"GGC",0b11101001},{"GGG",0b11101010},{"GGT",0b11101011},
	{"GTA",0b11101100},{"GTC",0b11101101},{"GTG",0b11101110},{"GTT",0b11101111},

	{"TAA",0b11110000},{"TAC",0b11110001},{"TAG",0b11110010},{"TAT",0b11110011},
	{"TCA",0b11110100},{"TCC",0b11110101},{"TCG",0b11110110},{"TCT",0b11110111},
	{"TGA",0b11111000},{"TGC",0b11111001},{"TGG",0b11111010},{"TGT",0b11111011},
	{"TTA",0b11111100},{"TTC",0b11111101},{"TTG",0b11111110},{"TTT",0b11111111},

	};




static const string int2string1[]={"A","C","G","T"};

static const string int2string2[]={
	"AA","AC","AG","AT",
	"CA","CC","CG","CT",
	"GA","GC","GG","GT",
	"TA","TC","TG","TT"
	};

static const string int2string3[]={
	"AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
	"CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
	"GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
	"TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"
	};


uint8_t char2Byte[]={
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,1,0,0,
	0,2,0,0,0,0,0,0,0,0,
	0,0,0,0,3
	};

uint8_t char2Byte1[]={
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,1<<2,0,0,
	0,2<<2,0,0,0,0,0,0,0,0,
	0,0,0,0,3<<2
	};

uint8_t char2Byte2[]={
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,1<<4,0,0,
	0,2<<4,0,0,0,0,0,0,0,0,
	0,0,0,0,3<<4
	};





void printUC(uint8_t a){
	for (int i = 0; i < 8; i++) {
		printf("%d", !!((a << i) & 0x80));
	}
	printf("\n");
}




binSeq::binSeq(const string& str){
	//~ isNumber=false;
	uint8_t mod(str.size()%3);
	size_t i(0);
	vect.reserve((str.size()/3)+mod);
	for (; i<str.size()-mod; i+=3){
		vect.push_back(char2int[functionwihtoutcoll(str[i],str[i+1],str[i+2])]);
		//~ vect.push_back(0x11000000 | char2Byte[str[i]] | char2Byte1[str[i+1]] | char2Byte2[str[i+2]]);
	}

	if(mod!=0){
		if(mod==1){
			vect.push_back(char2int1[str[i]]);
		}else{
			vect.push_back(char2int2[functionwihtoutcoll(str[i-1],str[i])]);
		}
	}
}



binSeq::binSeq(const binSeq& bs){
	//~ isNumber=bs.isNumber;
	vect=bs.vect;
}



string binSeq::str(){
	string res;
	res.reserve(vect.size());
	for(size_t i(0);i<vect.size();++i){

		switch (vect[i]>>6) {
			case 3:
				//~ res+=(int2string3[c]);
				res.append(int2string3[vect[i]&0b00111111]);
				break;
			case 2:
				res.append(int2string2[vect[i]&0b00001111]);
				//~ res+=(int2string2[c]);
				break;
			default:
				//~ res+=(int2string1[c]);
				res.append(int2string1[vect[i]&0b00000011]);
				break;
		}
	}
	return res;
}

void binSeq::str2(string& res){
	res.clear();
	res.reserve(vect.size());
	for(size_t i(0);i<vect.size();++i){

		switch (vect[i]>>6) {
			case 3:
				//~ res+=(int2string3[c]);
				res.append(int2string3[vect[i]&0b00111111]);
				break;
			case 2:
				res.append(int2string2[vect[i]&0b00001111]);
				//~ res+=(int2string2[c]);
				break;
			default:
				//~ res+=(int2string1[c]);
				res.append(int2string1[vect[i]&0b00000011]);
				break;
		}
	}
}


binSeq binSeq::sub(uint8_t begin){
	binSeq res;
	res.vect.reserve(vect.size());
	uint8_t count(0);
	bool go(true);
	size_t i(0);
	for(; go; ++i){
		uint8_t c(vect[i]);
		uint8_t n(c>>6);
		if(count+n<begin){
			count+=n;
		}else{
			go=false;
			if(count+n!=begin){
				if(count+n-begin==1){
					res.vect.push_back((1<<6)| c&0b00000011);
				}else{
					res.vect.push_back((2<<6)| c&0b00001111);
				}
			}
		}
	}
	res.vect.insert(res.vect.end(), vect.begin()+i, vect.end());
	return res;
}



//~ binSeq binSeq::sub(size_t begin, size_t size){
	//~ binSeq res;
	//~ res.vect.reserve(vect.size());
	//~ size_t i(0);
	//~ size_t countBegin(0);
	//~ bool go(false);
	//~ for (size_t count(0); count<size;) {
		//~ uint8_t ch(vect[i]);
		//~ uint8_t mod(ch/(1<<6));
//~
		//~ if(!go){
			//~ if(countBegin+mod<begin){
				//~ countBegin+=mod;
			//~ }else{
				//~ if(countBegin+mod==begin){
					//~ go=true;
				//~ }else{
					//~ go=true;
					//~ if(countBegin+mod-begin==1){
						//~ res.vect.push_back((1<<6)+ch%(1<<2));
						//~ ++count;
					//~ }else{
						//~ res.vect.push_back((2<<6)+ch%(1<<4));
						//~ count+=2;
					//~ }
				//~ }
			//~ }
		//~ }else{
			//~ if(count+mod<=size){
				//~ res.vect.push_back(ch);
				//~ count+=mod;
			//~ }else{
				//~ ch<<=2;
				//~ if(size-count==1){
					//~ res.vect.push_back((1<<6)+ch/(1<<6));
					//~ ++count;
				//~ }else{
					//~ res.vect.push_back((2<<6)+ch/(1<<4));
					//~ count+=2;
				//~ }
			//~ }
		//~ }
//~
		//~ ++i;
	//~ }
	//~ return res;
//~ }


binSeq binSeq::getBegin(uint8_t size){
	binSeq res;
	size_t i(0);
	for(uint8_t c(0); c<size;++i){
		uint8_t ch(vect[i]);
		uint8_t mod(ch/(1<<6));
		int8_t n(c+mod-size);

		switch (n){
			case 1:
				ch&=0b00111100;
				res.vect.push_back((2<<6)|(ch>>2));
				c+=2;
				break;
			case 2:
				ch&=0b00110000;
				res.vect.push_back((1<<6)|(ch>>4));
				c++;
				break;
			default:
				res.vect.push_back(ch);
				c+=mod;
				break;
		}
	}
	return res;
}

uint64_t binSeq::getBeginInt(uint8_t size){
	uint64_t res(0);
	size_t i(0);
	for(uint8_t c(0); c<size;++i){
		uint8_t ch(vect[i]);
		uint8_t mod(ch>>6);
		int8_t n(c+mod-size);

		switch (n){
			case 1:
				ch&=0b00111100;
				//~ res.vect.push_back((2<<6)|(ch>>2));
				res<<=4;
				res+=(ch>>2);
				c+=2;
				break;
			case 2:
				ch&=0b00110000;
				//~ res.vect.push_back((1<<6)|(ch>>4));
				res<<=2;
				res+=(ch>>4);
				c++;
				break;
			default:
				//~ res.vect.push_back(ch);
				res<<=(2*mod);
				res+=(ch&0b00111111);
				c+=mod;
				break;
		}
	}
	return res;
}

uint64_t binSeq::getBeginRcInt(uint8_t size){
	uint64_t res(0),inter;
	size_t i(0);
	for(uint8_t c(0); c<size;++i){
		uint8_t ch(rc[vect[i]]);
		uint8_t mod(ch>>6);
		int8_t n(c+mod-size);

		switch (n){
			case 1:
				inter=((ch&0b00001111));
				inter<<=(2*c);
				res+=inter;
				c+=2;
				break;
			case 2:
				inter=((ch&0b00000011));
				inter<<=(2*c);
				res+=inter;
				c++;
				break;
			default:
				inter=((ch&0b00111111));
				inter<<=(2*c);
				res+=inter;
				c+=mod;
				break;
		}
	}
	return res;
}


binSeq binSeq::getEnd(uint8_t size){
	binSeq res;
	size_t i(vect.size()-1);
	for(uint8_t c(0); c<size;--i){
		uint8_t ch(vect[i]);
		uint8_t mod(ch>>6);
		uint8_t n(c+mod-size);
		switch (n){
			case 1:
				res.vect.push_back((2<<6)|(ch&0b00001111));
				c+=2;
				break;
			case 2:
				res.vect.push_back((1<<6)|(ch&0b00000011));
				c++;
				break;
			default:
				res.vect.push_back(ch);
				c+=mod;
				break;
		}
	}
	::reverse(res.vect.begin(),res.vect.end());
	return res;
}

uint64_t binSeq::getEndInt(uint8_t size){
	uint64_t res(0),inter;
	size_t i(vect.size()-1);
	for(uint8_t c(0); c<size;--i){
		uint8_t ch(vect[i]);
		uint8_t mod(ch>>6);
		uint8_t n(c+mod-size);
		switch (n){
			case 1:
				//~ res.vect.push_back((2<<6)|(ch&0b00001111));
				inter=((ch&0b00001111));
				inter<<=(2*c);
				res+=inter;
				c+=2;
				break;
			case 2:
				//~ res.vect.push_back((1<<6)|(ch&0b00000011));
				inter=((ch&0b00000011));
				inter<<=(2*c);
				res+=inter;
				c++;
				break;
			default:
				//~ res.vect.push_back(ch);
				inter=((ch&0b00111111));
				inter<<=(2*c);
				res+=inter;
				c+=mod;
				break;
		}
		//~ printUC(c);
		//~ cout<<c<<" "<<res<<endl;
	}
	//~ ::reverse(res.vect.begin(),res.vect.end());
	//~ cin.get();
	return res;
}

uint64_t binSeq::getEndRcInt(uint8_t size){
	uint64_t res(0);
	size_t i(vect.size()-1);
	for(uint8_t c(0); c<size;--i){
		uint8_t ch(rc[vect[i]]);
		uint8_t mod(ch>>6);
		uint8_t n(c+mod-size);
		switch (n){
			case 1:
				ch&=0b00111100;
				res<<=4;
				res+=(ch>>2);
				c+=2;
				break;
			case 2:
				ch&=0b00110000;
				res<<=2;
				res+=(ch>>4);
				++c;
				break;
			default:
				res<<=(2*mod);
				res+=(ch&0b00111111);
				c+=mod;
				break;
		}
		//~ printUC(c);
		//~ cout<<c<<" "<<res<<endl;
	}
	//~ ::reverse(res.vect.begin(),res.vect.end());
	//~ cin.get();
	return res;
}



void binSeq::add(const binSeq& bs){
	vect.insert(vect.end(), bs.vect.begin(), bs.vect.end());
}


void binSeq::reverse(){
	uint32_t s(vect.size());
	uint32_t i(0);
	for(; i<s/2 ;++i){
		uint8_t inter (rc[vect[s-1-i]]);
		vect[s-1-i]=rc[vect[i]];
		vect[i]=inter;
	}
	if(s%2==1){
		vect[i]=rc[vect[i]];
	}
}


binSeq::binSeq(){
	//~ isNumber=false;
}


size_t binSeq::size(){
	return vect.size();
}


binSeq::binSeq(uint32_t n){
	//~ isNumber=true;
	while(n!=0){
		vect.push_back((uint8_t)n%(1<<8));
		n>>=8;
	}
}


uint32_t binSeq::getNumber(){
	uint32_t res(0);
	for(int i(vect.size()-1); i>-1; --i){
		res<<=8;
		res+=vect[i];
	}
	return res;
}


uint64_t binSeq::getInt(){
	uint64_t res(0);
	for(size_t i(0); i<vect.size(); ++i){
		uint8_t c(vect[i]);
		//~ uint8_t mod(c/(1<<6));
		res<<=(2*(c>>6));
		res+=c&0b00111111;
	}
	return res;
}


void binSeq::clear(){
	vect.clear();
}


binSeq binSeq::getReverse(){
	binSeq res;
	res.vect.reserve(vect.size());
	for(int i((int)vect.size()-1); i>-1; --i){
		res.vect.push_back(rc[vect[i]]);
	}
	return res;
}







































