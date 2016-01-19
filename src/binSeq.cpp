#include "binSeq.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <cmath>


uint min(uint a, uint b){return (a<b ? a : b);}


string RCnadine(const string& str){
	string res(str);
	uint n = str.size();
	for(int i(n-1), j(0); i > -1; i--, j++){
		unsigned char c = str[i];
		// unsigned char d = (c >> 4)&7;
		c ^= 4;
		if ((c&3) != 3)
			c ^= 17;
		res[j] = c;
		}
	return res;
}



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


static const string int2string[][64]={
	{},
	{"A","C","G","T"},
	{
		"AA","AC","AG","AT",
		"CA","CC","CG","CT",
		"GA","GC","GG","GT",
		"TA","TC","TG","TT"
	},
	{
		"AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
		"CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
		"GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
		"TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"
	}
};


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


uint functionwihtoutcoll(uint8_t x,uint8_t y,uint8_t z){
	return( (x&0b00000110)|((y&0b00000110)<<2)|((z&0b00000110)<<4) );
}


uint functionwihtoutcoll(uint8_t x,uint8_t y){
	return(x|(y<<4));
}


void initBinSeq(uint k){
	kmerBinseq=k-1;

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
	printf(" ");
}


binSeq::binSeq(const string& str){
	isInt=false;
	uint mod(str.size()%3);
	uint i(0);
	vect.reserve((str.size()/3)+1);
	for (; i<str.size()-mod; i+=3){
		vect.push_back(char2int[functionwihtoutcoll(str[i],str[i+1],str[i+2])]);
		//~ vect.push_back(0x11000000 | char2Byte[str[i]] | char2Byte1[str[i+1]] | char2Byte2[str[i+2]]);
	}
	if(mod!=0){
		if(mod==1){
			vect.push_back(char2int1[(uint8_t)(str[i])]);
		}else{
			vect.push_back(char2int2[functionwihtoutcoll(str[i-1],str[i])]);
		}
	}
    //std::cout << "vect size " << vect.size() << " capacity " << vect.capacity() << std::endl; // turns out it's not exactly as efficient as it should (size != capacity)
}


binSeq::binSeq(const binSeq& bs){
	isInt=false;
	vect=bs.vect;
}


string binSeq::str(){
	string res;
	res.reserve(3*vect.size());
	for(uint i(0);i<vect.size();++i){
		res.append(int2string[vect[i]>>6][vect[i]&0b00111111]);
	}
	// for(uint i(0);i<res.size();++i){
	// 	if(res[i]!='A' and res[i]!='C' and res[i]!='G' and res[i]!='T'){
	// 		cout<<"lol"<<endl;
	// 		exit(0);
	// 	}
	// }
	return res;
}


// void binSeq::str2(string& res){
// 	res.clear();
// 	res.reserve(3*vect.size());
// 	for(uint i(0);i<vect.size();++i){
// 		res.append(int2string[vect[i]>>6][vect[i]&0b00111111]);
// 	}
// }


// binSeq binSeq::sub(){
// 	uint begin(30);
// 	binSeq res;
// 	res.vect.reserve(vect.size());
// 	uint count(0);
// 	bool go(true);
// 	uint i(0);
// 	for(; go; ++i){
// 		uint c(vect[i]);
// 		uint n(c>>6);
// 		if(count+n<begin){
// 			count+=n;
// 		}else{
// 			go=false;
// 			if(count+n!=begin){
// 				if(count+n-begin==1){
// 					res.vect.push_back((1<<6)| c&0b00000011);
// 				}else{
// 					res.vect.push_back((2<<6)| c&0b00001111);
// 				}
// 			}
// 		}
// 	}
// 	res.vect.insert(res.vect.end(), vect.begin()+i, vect.end());
// 	return res;
// }



binSeq binSeq::sub(){
	binSeq res;
	res.vect.reserve(vect.size());
	uint count(0);
	bool go(true);
	uint i(0);
	for(; go; ++i){
		uint c(vect[i]);
		uint n(c>>6);
		if(count+n<kmerBinseq){
			count+=n;
		}else{
			go=false;
			if(count+n!=kmerBinseq){
				if(count+n-kmerBinseq==1){
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
//
// //
// binSeq binSeq::getBegin(uint size){
// 	binSeq res;
// 	uint i(0);
// 	for(uint c(0); c<size;++i){
// 		uint ch(vect[i]);
// 		uint mod(ch>>6);
// 		uint n(c+mod-size);
//
// 		switch (n){
// 			case 1:
// 				ch&=0b00111100;
// 				res.vect.push_back((2<<6)|(ch>>2));
// 				c+=2;
// 				break;
// 			case 2:
// 				ch&=0b00110000;
// 				res.vect.push_back((1<<6)|(ch>>4));
// 				c++;
// 				break;
// 			default:
// 				res.vect.push_back(ch);
// 				c+=mod;
// 				break;
// 		}
// 	}
// 	return res;
// }


kmer binSeq::getBeginInt(){
	kmer res(0);
	uint i(0);
	for(uint c(0); c<kmerBinseq;++i){
		uint ch(vect[i]);
		uint mod(ch>>6);
		int n(c+mod-kmerBinseq);
		if(n<=0){
			res<<=(2*mod);
			res+=(ch&0b00111111);
			c+=mod;
		}else{
			if(n==2){
				ch&=0b00110000;
				res<<=2;
				res+=(ch>>4);
				return res;
			}else{
				ch&=0b00111100;
				res<<=(2*mod-2);
				res+=(ch>>2);
				return res;
			}
		}
	}
	return res;
}

//
// kmer binSeq::getBeginRcInt(uint size){
// 	kmer res(0),inter;
// 	uint i(0);
// 	for(uint c(0); c<size;++i){
// 		uint ch(rc[vect[i]]);
// 		uint mod(ch>>6);
// 		int n(c+mod-size);
// 		if(n<=0){
// 			inter=(ch&0b00111111);
// 				inter<<=(2*c);
// 				res+=inter;
// 				c+=mod;
// 		}else{
// 			if(n==1){
// 				if(mod==3){
// 					inter=(ch&0b00001111);
// 				}else{
// 					inter=(ch&0b0000011);
// 				}
// 				inter<<=(2*c);
// 				res+=inter;
// 				return res;
// 			}else{
// 				inter=((ch&0b00000011));
// 				inter<<=(2*c);
// 				res+=inter;
// 				return res;
// 			}
// 		}
// 	}
// 	return res;
// }
//
// //
// binSeq binSeq::getEnd(uint size){
// 	binSeq res;
// 	uint i(vect.size()-1);
// 	for(uint c(0); c<size;--i){
// 		uint ch(vect[i]);
// 		uint mod(ch>>6);
// 		int n(c+mod-size);
// 		switch (n){
// 			case 1:
// 				res.vect.push_back((2<<6)|(ch&0b00001111));
// 				c+=2;
// 				break;
// 			case 2:
// 				res.vect.push_back((1<<6)|(ch&0b00000011));
// 				c++;
// 				break;
// 			default:
// 				res.vect.push_back(ch);
// 				c+=mod;
// 				break;
// 		}
// 	}
// 	::reverse(res.vect.begin(),res.vect.end());
// 	return res;
// }
//

kmer binSeq::getEndInt(){
	kmer res(0),inter(0);
	uint i(vect.size()-1);
	for(uint c(0); c<kmerBinseq;--i){
		uint ch(vect[i]);
		uint mod(ch>>6);
		int n(c+mod-kmerBinseq);
		if(n<=0){
			inter=(ch&0b00111111);
			inter<<=(2*c);
			res+=inter;
			c+=mod;
		}else{
			if(n==1){
				if(mod==2){
					inter=(ch&0b00000011);
				}else{
					inter=(ch&0b00001111);
				}
				inter<<=(2*c);
				res+=inter;
				//~ c+=2;
				return res;
			}else{
				inter=(ch&0b00000011);
				inter<<=(2*c);
				res+=inter;
				//~ c++;
				return res;
			}
		}
	}
	return res;
}

//
// kmer binSeq::getEndRcInt(uint size){
// 	kmer res(0);
// 	uint i(vect.size()-1);
// 	for(uint c(0); c<size;--i){
// 		uint ch(rc[vect[i]]);
// 		uint mod(ch>>6);
// 		int n(c+mod-size);
// 		if(n<=0){
// 			res<<=(2*mod);
// 			res+=(ch&0b00111111);
// 			c+=mod;
// 		}else{
// 			if(n==1){
// 				if(mod==3){
// 					ch&=0b00111100;
// 					res<<=4;
// 					res+=(ch>>2);
// 					return res;
// 				}else{
// 					ch&=0b00110000;
// 					res<<=2;
// 					res+=(ch>>4);
// 					return res;
// 				}
// 			}else{
// 				ch&=0b00110000;
// 				res<<=2;
// 				res+=(ch>>4);
// 				return res;
// 			}
// 		}
// 	}
// 	return res;
// }
//

void binSeq::add(const binSeq& bs){vect.insert(vect.end(), bs.vect.begin(), bs.vect.end());}


binSeq::binSeq(){vect.reserve(kmerBinseq/3+1);isInt=false;}


uint binSeq::size(){return vect.size();}


binSeq::binSeq(uint n){
	isInt=true;
	while(n!=0){
		vect.push_back((uint8_t)n%(1<<8));
		n>>=8;
	}
}


uint binSeq::getNumber(){
	uint res(0);
	for(int i(vect.size()-1); i>-1; --i){
		res<<=8;
		res+=vect[i];
	}
	return res;
}

//
// kmer binSeq::getInt(){
// 	kmer res(0);
// 	for(uint i(0); i<vect.size(); ++i){
// 		uint c(vect[i]);
// 		//~ uint mod(c/(1<<6));
// 		res<<=(2*(c>>6));
// 		res+=c&0b00111111;
// 	}
// 	return res;
// }


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


void binSeq::reverse(){
	//~ cout<<str()<<endl;
	// string lol(RCnadine(str()));
	uint s(vect.size());
	uint i(0);
	for(; i<s/2 ;++i){
		uint inter (rc[vect[s-1-i]]);
		vect[s-1-i]=rc[vect[i]];
		vect[i]=inter;
	}
	if(s%2==1){
		vect[i]=rc[vect[i]];
	}
}
