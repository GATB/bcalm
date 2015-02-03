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




const uint8_t rc[]={

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




const string int2string1[]={"A","C","G","T"};

const string int2string2[]={
	"AA","AC","AG","AT",
	"CA","CC","CG","CT",
	"GA","GC","GG","GT",
	"TA","TC","TG","TT"
	};

const string int2string3[]={
	"AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
	"CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
	"GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
	"TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"
	};

//I'm so sorry...
uint8_t char2Byte[85]={
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

	uint8_t char2Byte1[85]={
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

	uint8_t char2Byte2[85]={
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




//~ const unsigned char first



char revcompnadine (char s) {
	if (s == 'A') return 'T';
	else if (s == 'C') return 'G';
	else if (s == 'G') return 'C';
	else if (s == 'T') return 'A';
	else if (s == 'a') return 't';
	else if (s == 'c') return 'g';
	else if (s == 'g') return 'c';
	else if (s == 't') return 'a';
	return 'X';
}


//string revcomp (const string &s) {
//	string rc;
//	for (int i = (int)s.length() - 1; i >= 0; i--) rc += revcomp(s[i]);
//	return rc;
//}


string reversecomplementnadine (const string &s){
	string rc;
	for (int i = (int)s.length() - 1; i >= 0; i--) rc += revcompnadine(s[i]);
	return rc;
}






void printUC(uint8_t a){
	for (int i = 0; i < 8; i++) {
		printf("%d", !!((a << i) & 0x80));
	}
	printf("\n");

}



uint8_t char2int( uint8_t c){
	switch (c) {
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
	}
	//~ if(c<'G'){
		//~ if(c=='A'){
			//~ return 0;
		//~ }
		//~ return 1;
	//~ }
	//~ if(c=='G'){
		//~ return 2;
	//~ }
	//~ return 3;
	//~ return 3;
}

uint8_t char2int1( uint8_t c){
	switch (c) {
		case 'A':
			return 0;
		case 'C':
			return 1<<2;
		case 'G':
			return 2<<2;
		default:
			return 3<<2;
	}
	//~ if(c<'G'){
		//~ if(c=='A'){
			//~ return 0;
		//~ }
		//~ return 1<<2;
	//~ }
	//~ if(c=='G'){
		//~ return 2<<2;
	//~ }
	//~ return 3<<2;

	//~ return 3<<2;
}

uint8_t char2int2( uint8_t c){
	switch (c) {
		case 'A':
			return 0;
		case 'C':
			return 1<<4;
		case 'G':
			return 2<<4;
		default:
			return 3<<4;
	}
	//~ if(c<'G'){
		//~ if(c=='A'){
			//~ return 0;
		//~ }
		//~ return 1<<4;
	//~ }
	//~ if(c=='G'){
		//~ return 2<<4;
	//~ }
	//~ return 3<<4;
	//~ return 3<<4;
}



uint8_t int2char(uint8_t c){
	switch (c) {
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		default:
		//~ case 3:
			return 'T';
	}
	//~ return 'T';
//~
	//~ if(c<2){
		//~ if(c==0){
			//~ return 'A';
		//~ }
		//~ return 'C';
	//~ }
	//~ if(c==2){
		//~ return 'G';
	//~ }
	//~ return 'T';
}




binSeq::binSeq(const string& str){
	isNumber=false;
	uint8_t mod(str.size()%3);
	uint8_t c;
	size_t i(0);
	size_t j(0);
	vect.reserve(str.size()/3+mod);
	//~ vect.resize(str.size()/3+mod);
	for (; i<str.size()-mod; i+=3){
		//~ vect.push_back(string2Byte.at(str.substr(i,3)));

		//~ c=12;
		//~ c+=char2int(str[i]);
		//~ c<<=2;
		//~ c+=char2int(str[i+1]);
		//~ c<<=2;
		//~ c+=char2int(str[i+2]);
		c=(3<<6)+char2Byte2[str[i]]+char2Byte1[str[i+1]]+char2Byte[str[i+2]];

		vect.push_back(c);
		//~ vect[j++]=c;
		//~ if(c!=vect[vect.size()-1]){
			//~ printUC(vect[vect.size()-1]);
			//~ printUC(c);
			//~ cin.get();
		//~ }
	}

	//~ vect.push_back(string2Byte.at(str.substr(i)));

	//~ switch (mod) {
		//~ case 0:
			//~ break;
		//~ case 1:
			//c=1<<6;
			//c+=char2int(str[str.size()-1]);
//~
			//~ break;
		//~ case 2:
			//c=1<<5;
			//c+=char2int(str[str.size()-2]);
			//c<<=2;
			//~ c=1<<5+char2int1(str[str.size()-2])+char2int(str[str.size()-1]);
			//c+=char2int(str[str.size()-1]);
			//~ vect.push_back(c);
			//~ break;
	//~ }
	if(mod!=0){
		if(mod==1){
			c=(1<<6)+char2Byte[str[str.size()-1]];
			//~ vect[j]=c;
		}else{
			c=(1<<7)+char2Byte1[str[str.size()-2]]+char2Byte[str[str.size()-1]];
			//~ vect[j]=c;
		}
		vect.push_back(c);
		//~ vect[j]=c;
	}

}



binSeq::binSeq(const binSeq& bs){
	isNumber=bs.isNumber;
	vect=bs.vect;
}



string binSeq::str(){
	string res;
	res.reserve(vect.size());
//	cout<<vect.size()<<endl;
	for(size_t i(0);i<vect.size();++i){
		uint8_t c(vect[i]);
		uint8_t mod(c/(1<<6));
//		printUC(c);
		c%=1<<6;
//		cout<<"mod : ";
//		printUC(mod);

		switch (mod) {
			case 1:
				//~ res.push_back(int2char(c/(1<<2)));
				res+=(int2string1[c]);
				break;
			case 2:
				//~ res.push_back(int2char(c/(1<<4)));
				//~ c<<=2;
				//~ res.push_back(int2char((c/(1<<4))%4));
				res+=(int2string2[c]);
				break;

			default:
				//~ res.push_back(int2char(c/(1<<6)));
				//~ c<<=2;
				//~ res.push_back(int2char(c/(1<<6)));
				//~ c<<=2;
				//~ res.push_back(int2char(c/(1<<6)));
				res+=(int2string3[c]);
				break;
		}
	}
	return res;
}


binSeq binSeq::sub(size_t begin){
	binSeq res;
	res.vect.reserve(vect.size());
	size_t count(0);
	bool go(true);
	size_t i(0);
	for(; i<vect.size() and go; ++i){
		uint8_t c(vect[i]);
		uint8_t n(c/(1<<6));
		if(count+n<begin){
			count+=n;
		}else{
			go=false;
			if(count+n==begin){
			}else{
				//~ unsigned char toGet((unsigned char)());
				if(count+n-begin==1){
					res.vect.push_back((1<<6)+c%(1<<2));
				}else{
					res.vect.push_back((2<<6)+c%(1<<4));
				}
			}
		}
	}
	res.vect.insert(res.vect.end(), vect.begin()+i, vect.end());

	return res;
}


binSeq binSeq::sub(size_t begin, size_t size){
	binSeq res;
	res.vect.reserve(vect.size());
	size_t i(0);
	size_t countBegin(0);
	bool go(false);
	for (size_t count(0); count<size;) {
		uint8_t ch(vect[i]);
		uint8_t mod(ch/(1<<6));

		if(!go){
			if(countBegin+mod<begin){
				countBegin+=mod;
			}else{
				if(countBegin+mod==begin){
					go=true;
				}else{
					go=true;
					//~ unsigned char toGet((unsigned char)());
					if(countBegin+mod-begin==1){
						res.vect.push_back((1<<6)+ch%(1<<2));
						++count;
					}else{
						res.vect.push_back((2<<6)+ch%(1<<4));
						count+=2;
					}
				}
			}
		}else{
			if(count+mod<=size){
				res.vect.push_back(ch);
				count+=mod;
//				cout<<1<<endl;
			}else{
//				cout<<2<<endl;
				//~ unsigned char toGet((unsigned char)(size-count));
//				printUC(toGet);
				ch<<=2;
				if(size-count==1){
					res.vect.push_back((1<<6)+ch/(1<<6));
					++count;
				}else{
					res.vect.push_back((2<<6)+ch/(1<<4));
					count+=2;
				}
				//~ res.vect.push_back((unsigned char)(toGet<<6)+ch/(1<<(8-2*(toGet))));
//				printUC((unsigned char)(toGet<<6)+ch%(1<<(6-2*(3-toGet))));
				//~ count+=toGet;
			}
		}

		++i;
	}

	return res;
}


/*
uint64_t binSeq::getBegin(size_t size){
	uint64_t res(0);
	size_t i(0);
	for(size_t c(0); c<size;){
		unsigned char ch(vect[i]);
		unsigned char mod(ch/(1<<6));
		int64_t n(c+mod-size);
		if(n<=0){
			res<<=(2*mod);
			res+=ch%(1<<6);
			i++;
			c+=mod;
		}else{
			if(n==2){
				res<<=2;
				ch<<=2;
				res+=ch/(1<<6);
				c++;
				i++;
			}else{
				res<<=4;
				ch<<=2;
				res+=ch/(1<<4);
				c+=2;
				i++;
			}
		}
	}
	return res;
}
*/


binSeq binSeq::getBegin(size_t size){
	binSeq res;
	size_t i(0);
	for(size_t c(0); c<size;){
		uint8_t ch(vect[i]);
		uint8_t mod(ch/(1<<6));
		int64_t n(c+mod-size);
		if(n<=0){
//			res<<=(2*mod);
//			res+=ch%(1<<6);
			res.vect.push_back(ch);
			i++;
			c+=mod;
		}else{
			if(n==2){
//				res<<=2;
				ch<<=2;
//				res+=ch/(1<<6);
				res.vect.push_back((1<<6)+ch/(1<<6));
				c++;
				i++;
			}else{
//				res<<=4;
				ch<<=2;
//				res+=ch/(1<<4);
				res.vect.push_back((2<<6)+ch/(1<<4));
				c+=2;
				i++;
			}
		}
	}
	return res;
}


/*

uint64_t binSeq::getEnd(size_t size){
	uint64_t res(0);
	size_t i(vect.size()-1);
	for(size_t c(0); c<size;){
		unsigned char ch(vect[i]);
		unsigned char mod(ch/(1<<6));
		int64_t n(c+mod-size);
		if(n<=0){
			res+=(ch%(1<<6)<<(2*c));
			i--;
			c+=mod;
		}else{
			if(n==2){
				res+=(ch%(1<<2)<<(2*c));
				c++;
				i--;
			}else{
				res+=(ch%(1<<4)<<(2*c));
				c+=2;
				i--;
			}
		}
	}
	return res;
}
*/

binSeq binSeq::getEnd(size_t size){
	binSeq res;
	size_t i(vect.size()-1);
	for(size_t c(0); c<size;){
		uint8_t ch(vect[i]);
		uint8_t mod(ch/(1<<6));
		int64_t n(c+mod-size);
		if(n<=0){
//			res+=(ch%(1<<6)<<(2*c));
			res.vect.push_back(ch);
			i--;
			c+=mod;
		}else{
			if(n==2){
//				res+=(ch%(1<<2)<<(2*c));
				res.vect.push_back((1<<6)+ch%(1<<2));
				c++;
				i--;
			}else{
//				res+=(ch%(1<<4)<<(2*c));
				res.vect.push_back((2<<6)+ch%(1<<4));
				c+=2;
				i--;
			}
		}
	}
	::reverse(res.vect.begin(),res.vect.end());
	return res;
}



void binSeq::add(const binSeq& bs){
	vect.insert(vect.end(), bs.vect.begin(), bs.vect.end());
}

//~
//~ void binSeq::reverse(){//TODO IMPROVE IN PLACE
	//~ vector<uint8_t> V;
	//~ V.reserve(vect.size());
	//~ for(int i((int)vect.size()-1); i>-1; --i){
		//~ V.push_back(rc[vect[i]]);
	//~ }
	//~ vect=V;
//~ }

void binSeq::reverse(){
	uint32_t s(vect.size());
	string str1(str());
	uint32_t i(0);
	for(; i<s/2 ;++i){
		uint8_t inter (rc[vect[s-1-i]]);
		vect[s-1-i]=rc[vect[i]];
		vect[i]=inter;
	}
	if(s%2==1){
		vect[i]=rc[vect[i]];
	}
	//~ if(str()!=reversecomplementnadine(str1)){
		//~ cout<<str1<<" "<<str()<<endl;
	//~ }
}


binSeq::binSeq(){
	isNumber=false;
}


size_t binSeq::size(){
	return vect.size();
}


binSeq::binSeq(uint32_t n){
	isNumber=true;
	while(n!=0){
		vect.push_back((uint8_t)n%(1<<8));
		n>>=8;
	}
}


uint32_t binSeq::getNumber(){
	uint32_t res(0);
	for(int i(vect.size()-1); i>-1; --i){
		uint8_t c(vect[i]);
		res<<=8;
		res+=c;
	}
	return res;
}


uint64_t binSeq::getInt(){
	uint64_t res(0);
	for(size_t i(0); i<vect.size(); ++i){
		uint8_t c(vect[i]);
		uint8_t mod(c/(1<<6));
		res<<=(2*mod);
		res+=c%(1<<(2*mod));
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







































