#include <glue.hpp>
#include <assert.h>
#include <gatb/gatb_core.hpp>
#include <ograph.h>
#include <iostream>
#include <memory>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#ifndef OSX
#include <sys/sysinfo.h> // to determine system memory
#endif

//#define SPARSEHASH 
#ifdef SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;
#endif 

string key_of_interest = "somestringthatwillneveroccur";
//string key_of_interest = "GATGTGTTTGTCGTGTTCGCGATCGCCACGGCATCGGTTAGACTGCTAACC";
bool glueDebug = false;
bool weWantSpeed = false;
using namespace std;

uint8_t byteLookupTable[256][7] = {{'A', 'A', 'A', 'A', '0', '0', '0'}, {'A', 'A', 'A', 'C', '0', '0', '0'}, {'A', 'A', 'A', 'G', '0', '0', '0'}, {'A', 'A', 'A', 'T', '0', '0', '0'}, {'A', 'A', 'C', 'A', '0', '0', '0'}, {'A', 'A', 'C', 'C', '0', '0', '0'}, {'A', 'A', 'C', 'G', '0', '0', '0'}, {'A', 'A', 'C', 'T', '0', '0', '0'}, {'A', 'A', 'G', 'A', '0', '0', '0'}, {'A', 'A', 'G', 'C', '0', '0', '0'}, {'A', 'A', 'G', 'G', '0', '0', '0'}, {'A', 'A', 'G', 'T', '0', '0', '0'}, {'A', 'A', 'T', 'A', '0', '0', '0'}, {'A', 'A', 'T', 'C', '0', '0', '0'}, {'A', 'A', 'T', 'G', '0', '0', '0'}, {'A', 'A', 'T', 'T', '0', '0', '0'}, {'A', 'C', 'A', 'A', '0', '0', '2'}, {'A', 'C', 'A', 'C', '0', '0', '2'}, {'A', 'C', 'A', 'G', '0', '0', '2'}, {'A', 'C', 'A', 'T', '0', '0', '2'}, {'A', 'C', 'C', 'A', '0', '0', '2'}, {'A', 'C', 'C', 'C', '0', '0', '2'}, {'A', 'C', 'C', 'G', '0', '0', '2'}, {'A', 'C', 'C', 'T', '0', '0', '2'}, {'A', 'C', 'G', 'A', '0', '0', '2'}, {'A', 'C', 'G', 'C', '0', '0', '2'}, {'A', 'C', 'G', 'G', '0', '0', '2'}, {'A', 'C', 'G', 'T', '0', '0', '2'}, {'A', 'C', 'T', 'A', '0', '0', '2'}, {'A', 'C', 'T', 'C', '0', '0', '2'}, {'A', 'C', 'T', 'G', '0', '0', '2'}, {'A', 'C', 'T', 'T', '0', '0', '2'}, {'A', 'G', 'A', 'A', '0', '0', '4'}, {'A', 'G', 'A', 'C', '0', '0', '4'}, {'A', 'G', 'A', 'G', '0', '0', '4'}, {'A', 'G', 'A', 'T', '0', '0', '4'}, {'A', 'G', 'C', 'A', '0', '0', '4'}, {'A', 'G', 'C', 'C', '0', '0', '4'}, {'A', 'G', 'C', 'G', '0', '0', '4'}, {'A', 'G', 'C', 'T', '0', '0', '4'}, {'A', 'G', 'G', 'A', '0', '0', '4'}, {'A', 'G', 'G', 'C', '0', '0', '4'}, {'A', 'G', 'G', 'G', '0', '0', '4'}, {'A', 'G', 'G', 'T', '0', '0', '4'}, {'A', 'G', 'T', 'A', '0', '0', '4'}, {'A', 'G', 'T', 'C', '0', '0', '4'}, {'A', 'G', 'T', 'G', '0', '0', '4'}, {'A', 'G', 'T', 'T', '0', '0', '4'}, {'A', 'T', 'A', 'A', '0', '0', '6'}, {'A', 'T', 'A', 'C', '0', '0', '6'}, {'A', 'T', 'A', 'G', '0', '0', '6'}, {'A', 'T', 'A', 'T', '0', '0', '6'}, {'A', 'T', 'C', 'A', '0', '0', '6'}, {'A', 'T', 'C', 'C', '0', '0', '6'}, {'A', 'T', 'C', 'G', '0', '0', '6'}, {'A', 'T', 'C', 'T', '0', '0', '6'}, {'A', 'T', 'G', 'A', '0', '0', '6'}, {'A', 'T', 'G', 'C', '0', '0', '6'}, {'A', 'T', 'G', 'G', '0', '0', '6'}, {'A', 'T', 'G', 'T', '0', '0', '6'}, {'A', 'T', 'T', 'A', '0', '0', '6'}, {'A', 'T', 'T', 'C', '0', '0', '6'}, {'A', 'T', 'T', 'G', '0', '0', '6'}, {'A', 'T', 'T', 'T', '0', '0', '6'}, {'C', 'A', 'A', 'A', '0', '1', '0'}, {'C', 'A', 'A', 'C', '0', '1', '0'}, {'C', 'A', 'A', 'G', '0', '1', '0'}, {'C', 'A', 'A', 'T', '0', '1', '0'}, {'C', 'A', 'C', 'A', '0', '1', '0'}, {'C', 'A', 'C', 'C', '0', '1', '0'}, {'C', 'A', 'C', 'G', '0', '1', '0'}, {'C', 'A', 'C', 'T', '0', '1', '0'}, {'C', 'A', 'G', 'A', '0', '1', '0'}, {'C', 'A', 'G', 'C', '0', '1', '0'}, {'C', 'A', 'G', 'G', '0', '1', '0'}, {'C', 'A', 'G', 'T', '0', '1', '0'}, {'C', 'A', 'T', 'A', '0', '1', '0'}, {'C', 'A', 'T', 'C', '0', '1', '0'}, {'C', 'A', 'T', 'G', '0', '1', '0'}, {'C', 'A', 'T', 'T', '0', '1', '0'}, {'C', 'C', 'A', 'A', '0', '1', '2'}, {'C', 'C', 'A', 'C', '0', '1', '2'}, {'C', 'C', 'A', 'G', '0', '1', '2'}, {'C', 'C', 'A', 'T', '0', '1', '2'}, {'C', 'C', 'C', 'A', '0', '1', '2'}, {'C', 'C', 'C', 'C', '0', '1', '2'}, {'C', 'C', 'C', 'G', '0', '1', '2'}, {'C', 'C', 'C', 'T', '0', '1', '2'}, {'C', 'C', 'G', 'A', '0', '1', '2'}, {'C', 'C', 'G', 'C', '0', '1', '2'}, {'C', 'C', 'G', 'G', '0', '1', '2'}, {'C', 'C', 'G', 'T', '0', '1', '2'}, {'C', 'C', 'T', 'A', '0', '1', '2'}, {'C', 'C', 'T', 'C', '0', '1', '2'}, {'C', 'C', 'T', 'G', '0', '1', '2'}, {'C', 'C', 'T', 'T', '0', '1', '2'}, {'C', 'G', 'A', 'A', '0', '1', '4'}, {'C', 'G', 'A', 'C', '0', '1', '4'}, {'C', 'G', 'A', 'G', '0', '1', '4'}, {'C', 'G', 'A', 'T', '0', '1', '4'}, {'C', 'G', 'C', 'A', '0', '1', '4'}, {'C', 'G', 'C', 'C', '0', '1', '4'}, {'C', 'G', 'C', 'G', '0', '1', '4'}, {'C', 'G', 'C', 'T', '0', '1', '4'}, {'C', 'G', 'G', 'A', '0', '1', '4'}, {'C', 'G', 'G', 'C', '0', '1', '4'}, {'C', 'G', 'G', 'G', '0', '1', '4'}, {'C', 'G', 'G', 'T', '0', '1', '4'}, {'C', 'G', 'T', 'A', '0', '1', '4'}, {'C', 'G', 'T', 'C', '0', '1', '4'}, {'C', 'G', 'T', 'G', '0', '1', '4'}, {'C', 'G', 'T', 'T', '0', '1', '4'}, {'C', 'T', 'A', 'A', '0', '1', '6'}, {'C', 'T', 'A', 'C', '0', '1', '6'}, {'C', 'T', 'A', 'G', '0', '1', '6'}, {'C', 'T', 'A', 'T', '0', '1', '6'}, {'C', 'T', 'C', 'A', '0', '1', '6'}, {'C', 'T', 'C', 'C', '0', '1', '6'}, {'C', 'T', 'C', 'G', '0', '1', '6'}, {'C', 'T', 'C', 'T', '0', '1', '6'}, {'C', 'T', 'G', 'A', '0', '1', '6'}, {'C', 'T', 'G', 'C', '0', '1', '6'}, {'C', 'T', 'G', 'G', '0', '1', '6'}, {'C', 'T', 'G', 'T', '0', '1', '6'}, {'C', 'T', 'T', 'A', '0', '1', '6'}, {'C', 'T', 'T', 'C', '0', '1', '6'}, {'C', 'T', 'T', 'G', '0', '1', '6'}, {'C', 'T', 'T', 'T', '0', '1', '6'}, {'G', 'A', 'A', 'A', '1', '0', '0'}, {'G', 'A', 'A', 'C', '1', '0', '0'}, {'G', 'A', 'A', 'G', '1', '0', '0'}, {'G', 'A', 'A', 'T', '1', '0', '0'}, {'G', 'A', 'C', 'A', '1', '0', '0'}, {'G', 'A', 'C', 'C', '1', '0', '0'}, {'G', 'A', 'C', 'G', '1', '0', '0'}, {'G', 'A', 'C', 'T', '1', '0', '0'}, {'G', 'A', 'G', 'A', '1', '0', '0'}, {'G', 'A', 'G', 'C', '1', '0', '0'}, {'G', 'A', 'G', 'G', '1', '0', '0'}, {'G', 'A', 'G', 'T', '1', '0', '0'}, {'G', 'A', 'T', 'A', '1', '0', '0'}, {'G', 'A', 'T', 'C', '1', '0', '0'}, {'G', 'A', 'T', 'G', '1', '0', '0'}, {'G', 'A', 'T', 'T', '1', '0', '0'}, {'G', 'C', 'A', 'A', '1', '0', '2'}, {'G', 'C', 'A', 'C', '1', '0', '2'}, {'G', 'C', 'A', 'G', '1', '0', '2'}, {'G', 'C', 'A', 'T', '1', '0', '2'}, {'G', 'C', 'C', 'A', '1', '0', '2'}, {'G', 'C', 'C', 'C', '1', '0', '2'}, {'G', 'C', 'C', 'G', '1', '0', '2'}, {'G', 'C', 'C', 'T', '1', '0', '2'}, {'G', 'C', 'G', 'A', '1', '0', '2'}, {'G', 'C', 'G', 'C', '1', '0', '2'}, {'G', 'C', 'G', 'G', '1', '0', '2'}, {'G', 'C', 'G', 'T', '1', '0', '2'}, {'G', 'C', 'T', 'A', '1', '0', '2'}, {'G', 'C', 'T', 'C', '1', '0', '2'}, {'G', 'C', 'T', 'G', '1', '0', '2'}, {'G', 'C', 'T', 'T', '1', '0', '2'}, {'G', 'G', 'A', 'A', '1', '0', '4'}, {'G', 'G', 'A', 'C', '1', '0', '4'}, {'G', 'G', 'A', 'G', '1', '0', '4'}, {'G', 'G', 'A', 'T', '1', '0', '4'}, {'G', 'G', 'C', 'A', '1', '0', '4'}, {'G', 'G', 'C', 'C', '1', '0', '4'}, {'G', 'G', 'C', 'G', '1', '0', '4'}, {'G', 'G', 'C', 'T', '1', '0', '4'}, {'G', 'G', 'G', 'A', '1', '0', '4'}, {'G', 'G', 'G', 'C', '1', '0', '4'}, {'G', 'G', 'G', 'G', '1', '0', '4'}, {'G', 'G', 'G', 'T', '1', '0', '4'}, {'G', 'G', 'T', 'A', '1', '0', '4'}, {'G', 'G', 'T', 'C', '1', '0', '4'}, {'G', 'G', 'T', 'G', '1', '0', '4'}, {'G', 'G', 'T', 'T', '1', '0', '4'}, {'G', 'T', 'A', 'A', '1', '0', '6'}, {'G', 'T', 'A', 'C', '1', '0', '6'}, {'G', 'T', 'A', 'G', '1', '0', '6'}, {'G', 'T', 'A', 'T', '1', '0', '6'}, {'G', 'T', 'C', 'A', '1', '0', '6'}, {'G', 'T', 'C', 'C', '1', '0', '6'}, {'G', 'T', 'C', 'G', '1', '0', '6'}, {'G', 'T', 'C', 'T', '1', '0', '6'}, {'G', 'T', 'G', 'A', '1', '0', '6'}, {'G', 'T', 'G', 'C', '1', '0', '6'}, {'G', 'T', 'G', 'G', '1', '0', '6'}, {'G', 'T', 'G', 'T', '1', '0', '6'}, {'G', 'T', 'T', 'A', '1', '0', '6'}, {'G', 'T', 'T', 'C', '1', '0', '6'}, {'G', 'T', 'T', 'G', '1', '0', '6'}, {'G', 'T', 'T', 'T', '1', '0', '6'}, {'T', 'A', 'A', 'A', '1', '1', '0'}, {'T', 'A', 'A', 'C', '1', '1', '0'}, {'T', 'A', 'A', 'G', '1', '1', '0'}, {'T', 'A', 'A', 'T', '1', '1', '0'}, {'T', 'A', 'C', 'A', '1', '1', '0'}, {'T', 'A', 'C', 'C', '1', '1', '0'}, {'T', 'A', 'C', 'G', '1', '1', '0'}, {'T', 'A', 'C', 'T', '1', '1', '0'}, {'T', 'A', 'G', 'A', '1', '1', '0'}, {'T', 'A', 'G', 'C', '1', '1', '0'}, {'T', 'A', 'G', 'G', '1', '1', '0'}, {'T', 'A', 'G', 'T', '1', '1', '0'}, {'T', 'A', 'T', 'A', '1', '1', '0'}, {'T', 'A', 'T', 'C', '1', '1', '0'}, {'T', 'A', 'T', 'G', '1', '1', '0'}, {'T', 'A', 'T', 'T', '1', '1', '0'}, {'T', 'C', 'A', 'A', '1', '1', '2'}, {'T', 'C', 'A', 'C', '1', '1', '2'}, {'T', 'C', 'A', 'G', '1', '1', '2'}, {'T', 'C', 'A', 'T', '1', '1', '2'}, {'T', 'C', 'C', 'A', '1', '1', '2'}, {'T', 'C', 'C', 'C', '1', '1', '2'}, {'T', 'C', 'C', 'G', '1', '1', '2'}, {'T', 'C', 'C', 'T', '1', '1', '2'}, {'T', 'C', 'G', 'A', '1', '1', '2'}, {'T', 'C', 'G', 'C', '1', '1', '2'}, {'T', 'C', 'G', 'G', '1', '1', '2'}, {'T', 'C', 'G', 'T', '1', '1', '2'}, {'T', 'C', 'T', 'A', '1', '1', '2'}, {'T', 'C', 'T', 'C', '1', '1', '2'}, {'T', 'C', 'T', 'G', '1', '1', '2'}, {'T', 'C', 'T', 'T', '1', '1', '2'}, {'T', 'G', 'A', 'A', '1', '1', '4'}, {'T', 'G', 'A', 'C', '1', '1', '4'}, {'T', 'G', 'A', 'G', '1', '1', '4'}, {'T', 'G', 'A', 'T', '1', '1', '4'}, {'T', 'G', 'C', 'A', '1', '1', '4'}, {'T', 'G', 'C', 'C', '1', '1', '4'}, {'T', 'G', 'C', 'G', '1', '1', '4'}, {'T', 'G', 'C', 'T', '1', '1', '4'}, {'T', 'G', 'G', 'A', '1', '1', '4'}, {'T', 'G', 'G', 'C', '1', '1', '4'}, {'T', 'G', 'G', 'G', '1', '1', '4'}, {'T', 'G', 'G', 'T', '1', '1', '4'}, {'T', 'G', 'T', 'A', '1', '1', '4'}, {'T', 'G', 'T', 'C', '1', '1', '4'}, {'T', 'G', 'T', 'G', '1', '1', '4'}, {'T', 'G', 'T', 'T', '1', '1', '4'}, {'T', 'T', 'A', 'A', '1', '1', '6'}, {'T', 'T', 'A', 'C', '1', '1', '6'}, {'T', 'T', 'A', 'G', '1', '1', '6'}, {'T', 'T', 'A', 'T', '1', '1', '6'}, {'T', 'T', 'C', 'A', '1', '1', '6'}, {'T', 'T', 'C', 'C', '1', '1', '6'}, {'T', 'T', 'C', 'G', '1', '1', '6'}, {'T', 'T', 'C', 'T', '1', '1', '6'}, {'T', 'T', 'G', 'A', '1', '1', '6'}, {'T', 'T', 'G', 'C', '1', '1', '6'}, {'T', 'T', 'G', 'G', '1', '1', '6'}, {'T', 'T', 'G', 'T', '1', '1', '6'}, {'T', 'T', 'T', 'A', '1', '1', '6'}, {'T', 'T', 'T', 'C', '1', '1', '6'}, {'T', 'T', 'T', 'G', '1', '1', '6'}, {'T', 'T', 'T', 'T', '1', '1', '6'}};

uint8_t bits2byteTable[2][2][2][2][2][2][2][2] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};

/*uint8_t bits2byte (RawEntry bits, int start) {
	return  
		1 * bits[start + 7] +  
		2 * bits[start + 6] + 
		4 * bits[start + 5] + 
		8 * bits[start + 4] + 
		16 * bits[start + 3] + 
		32 * bits[start + 2] + 
		64 * bits[start + 1] + 
		128 * bits[start];
}*/




template<typename T>
string add_commas(T num) {
	string s, retval;
	ostringstream o;
	o << num;
	s = o.str();
	for (int i = 0; i < s.length(); i++) {
		retval.push_back(s.at(i));
		int j = s.length() - 1 - i;
		if (((j % 3) == 0) && (j != 0)) {
			retval.push_back(',');
		}
	}
	return retval;
}


 //print vector of T
template<class T>
ostream & operator << (ostream & out, const vector<T> & v) {
	out << v.size() << '\t';
	if (v.size() != 0) {
		out << int(v[0]);
		for (int i = 1; i < v.size(); i++) out << '\t' << int(v[i]);
	}
	return out;
}

char revcomp (char s) {
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

string revcomp (string &s) {
	string rc;
	for (int i = s.length() - 1; i >= 0; i--) rc += revcomp(s[i]);
	return rc;
}

//this short-circuits the computation as soon as the difference is found
string rcnorm (string &s) {
	string rc;
	for (int i = 0; i < s.length(); i++) {
		char c = revcomp(s[s.length() - 1 - i]);
		if (s[i] < c) {
			return s;
		} else if (s[i] > c) {
			return revcomp(s);
		}
	}
	return revcomp(s);
}

//string rcnorm ( string seq ) { return std::min (seq, reversecomplement(seq)); }

/*char bits2char (bool leftbit, bool rightbit) {
	if (leftbit == 0 && rightbit == 0) return 'A';
	else if (leftbit == 0 && rightbit == 1) return 'C';
	else if (leftbit == 1 && rightbit == 0) return 'G';
	else return 'T';
}
*/

/* nifty function to highlight a substring for printing to terminal */
string debug_highlight(string s, string motif)
{
    size_t pos = s.find(motif);
    if (pos == string::npos)
        pos = s.find(revcomp(motif));
    if (pos == string::npos)
    {
        return s;
    }

    return s.substr(0,pos) + "\033[1;31m" + s.substr(pos, motif.size()) + "\033[0m" + s.substr(pos + motif.size(), s.size() - pos - motif.size()); 
}


string tostring(GlueEntry e, string key) {
	ostringstream out;
	if (e.getLmark()) out << "_"; else out << " ";
	out << debug_highlight(e.getSeq(), key);
	if (e.getRmark()) out << "_"; else out << " ";
	return out.str();
}

GlueEntry::GlueEntry(RawEntry raw, size_t _kmerSize){
	kmerSize = _kmerSize;
#ifdef TWOBITGLUEHASH

	//cout << "GlueEntry constructor got raw = " << raw << ".\n";
	uint8_t * vals = byteLookupTable[raw[0]];
	//cout << vals[0] << vals[1] << vals[2] << vals[3] << vals[4] << vals[5] << vals[6] << "\t:" << seq << ":" << endl; 
	lmark = vals[4] - '0';
	rmark = vals[5] - '0';
	int leftoverbits = vals[6] - '0';

	int seqLen;
	if (leftoverbits == 0) {
		seqLen = raw.size() * 4 - 2;
	} else {
		seqLen = raw.size() * 4 - 2 - (4 - (leftoverbits / 2));
	}
	seq.resize(seqLen);

	int seqCounter = 0;
	if (raw.size() == 1) {
		if (leftoverbits == 0) {
			seq[seqCounter++] = vals[2];
			seq[seqCounter++] = vals[3];
			//seq.push_back(vals[2]); seq.push_back(vals[3]);
		} else if (leftoverbits == 6) {
			seq[seqCounter++] = vals[2];
			//seq.push_back(vals[2]);
		} 
	} else {
		seq[seqCounter++] = vals[2];
		seq[seqCounter++] = vals[3];
		for (int i = 1; i < raw.size() - 1; i++) {
			vals =  byteLookupTable[raw[i]];
			seq[seqCounter++] = vals[0];
			seq[seqCounter++] = vals[1];
			seq[seqCounter++] = vals[2];
			seq[seqCounter++] = vals[3];
		}

		vals =  byteLookupTable[raw[raw.size() - 1]];
		if (leftoverbits == 0) {
			seq[seqCounter++] = vals[0];
			seq[seqCounter++] = vals[1];
			seq[seqCounter++] = vals[2];
			seq[seqCounter++] = vals[3];
		} else if (leftoverbits == 6) {
			seq[seqCounter++] = vals[0];
			seq[seqCounter++] = vals[1];
			seq[seqCounter++] = vals[2];
		} else if (leftoverbits == 4) {
			seq[seqCounter++] = vals[0];
			seq[seqCounter++] = vals[1];
		} else if (leftoverbits == 2) {
			seq[seqCounter++] = vals[0];
		}

	}
	//cout << "Constructor put togegher: " << seq << ".\n";

#else
	seq = raw.substr(0, raw.size() - 2);
	lmark = raw.at(raw.length() - 2);
	rmark = raw.at(raw.length() - 1);
#endif
}




RawEntry GlueEntry::getRaw() {
#ifdef TWOBITGLUEHASH
	//cout << "Entering getRaw with: " << tostring(*this, "CCCCC") << ".\n";
	int rawLengthBits = 0;
	rawLengthBits += 2; //to represent the left and right mark
	rawLengthBits += 2; //to represent how many bits are to be read from the last byte: 00 = 0 bits, 01 = 2 bits, 10 = 4 bits, 11 = 6 bits
	rawLengthBits += (seq.length() * 2); //to represent the sequence
	int leftoverbits = rawLengthBits - ((rawLengthBits / 8) * 8); 
	int rawLengthBytes = (rawLengthBits + 7.99) / 8;
	vector<uint8_t> bits(rawLengthBytes * 8, 0);
	bits[0] = lmark;
	bits[1] = rmark;
	if (leftoverbits == 0) {
		bits[2] = 0; bits[3] = 0;
	} else if (leftoverbits == 2) {
		bits[2] = 0; bits[3] = 1;
	} else if (leftoverbits == 4) {
		bits[2] = 1; bits[3] = 0;
	} else if (leftoverbits == 6) {
		bits[2] = 1; bits[3] = 1;
	} else {
		cout << "GlueEntry::getRaw(): Unexpected leftover bits: " << leftoverbits << ".\n";
		exit(1);
	}
	for (int i = 0; i < seq.length(); i++) {
		if ((seq[i] == 'A') || (seq[i] == 'C')) {
			bits[4 + 2*i] = 0;
		} else {
			bits[4 + 2*i] = 1;
		}
		if ((seq[i] == 'A') || (seq[i] == 'G')) {
			bits[4 + 2*i + 1] = 0;
		} else {
			bits[4 + 2*i + 1] = 1;
		}
	}
	//encode bits in a 2bit compacted string.
	RawEntry raw;
	for (int i = 0; i < bits.size(); i += 8) {
		raw.push_back(bits2byteTable[bits[i]][bits[i+1]][bits[i+2]][bits[i+3]][bits[i+4]][bits[i+5]][bits[i+6]][bits[i+7]]);
	}
	//cout << "bits: " << bits << ".\n";
	return raw;

#else
	string raw;
	raw = seq;
	raw.push_back(lmark);
	raw.push_back(rmark);
#endif

	return raw;
}


string GlueStorage::dump() { 
	ostringstream o;
	for (auto it = glueMap.begin(); it != glueMap.end(); it++) {
		GlueEntry e;
		derefIt(it, e);
		string key = it->first;
		o << "\"" << key << "\"\t\"" << tostring(e, key) << "\"" << endl;
	}
	return o.str();
}

string GlueStorage::dump(string key, bool dumpkey) { 
	ostringstream o;
	GlueEntry e;
	if(find(key, e)) {
		if (dumpkey) o << "\"" << key << "\"\t";
		o << "\"" << tostring(e, key) << "\"";
	}
	return o.str();
}


bool GlueStorage::derefIt (GlueMap::const_iterator it, GlueEntry & e) {
	if (it == glueMap.end()) return false;
	e = GlueEntry(it->second, kmerSize); 
	return true;
}

bool GlueStorage::find (string key, GlueEntry & e) {
	findIt  = glueMap.find (key);
	return derefIt(findIt, e);
}

void GlueStorage::insertAtKey(string key, GlueEntry e) {
	glueMap[key] = e.getRaw();
}

void GlueStorage::insertAtIt(GlueMap::iterator it, GlueEntry e) {
	if (it == glueMap.end()) {
		cout << "GlueStorage::insertAtIt has it pointing to the end of glueMap...SHOULD not BE.\n";
		exit(1);
	}
	it->second = e.getRaw();
}

//glueResult contains the result of the glue
void Glue::glueSingleEntry(GlueEntry query, GlueEntry match, string key, GlueEntry & glueResult) {
	glueResult = GlueEntry();
	
	if (weWantSpeed) { //this disables the error check, meaning that if there is some bug in the program that causes query and match to be non-gluable, the program will not terminate and will continue
		if (match.getLmark()  && (query.getRkmer() == match.getLkmer()) && (query.getRkmer() == key)) {
			string gluedStr = query.getSeq() + match.getSeq().substr(kmerSize, match.getSeq().size() - kmerSize);
			glueResult = GlueEntry(gluedStr, query.getLmark(), match.getRmark(), kmerSize);
		} else {
			string gluedStr = match.getSeq() + query.getSeq().substr(kmerSize, query.getSeq().size() - kmerSize);
			glueResult = GlueEntry(gluedStr, match.getLmark(), query.getRmark(), kmerSize);
		}
	} else{
		if (match.getLmark()  && (query.getRkmer() == match.getLkmer()) && (query.getRkmer() == key)) {
			string gluedStr = query.getSeq() + match.getSeq().substr(kmerSize, match.getSeq().size() - kmerSize);
			glueResult = GlueEntry(gluedStr, query.getLmark(), match.getRmark(), kmerSize);
		} else if (query.getLmark()  && (match.getRkmer() == query.getLkmer()) && (match.getRkmer() == key)) {
			string gluedStr = match.getSeq() + query.getSeq().substr(kmerSize, query.getSeq().size() - kmerSize);
			glueResult = GlueEntry(gluedStr, match.getLmark(), query.getRmark(), kmerSize);
		} else {
			cout << "glueSingleEntry received non-gluable input:\n";
			cout << "  query: " << tostring(query, key) << endl;
			cout << "  match: " << tostring(match, key) << endl;
			exit(1);
		}
	}
	
}


//returns true if inserted
bool Glue::insert_aux(GlueEntry newEntry, string key_norm, GlueEntry & glueResult, bool onlyInsertIfFull) { 
	GlueEntry e;

	//if (key == key_of_interest)  cout << key << "\tinsert_aux\t" << tostring(newEntry, key)  << "\tprior\t" << glueStorage.dump(key, false) << "\t";
	if (glueDebug) cout << key_norm << "\tinsert_aux\t" << tostring(newEntry, key_norm)  << "\tprior\t" << glueStorage.dump(key_norm, false) << "\t";

	if (!glueStorage.find(key_norm, e)) {
		if (onlyInsertIfFull) return false;
		glueStorage.insertAtKey(key_norm, newEntry);
	} else {
		if (e.isEmpty()) {
			if (onlyInsertIfFull) return false;
			glueStorage.insertAfterFind(newEntry);
			//glueStorage.insertAtKey(key_norm, newEntry);
		} else {
			//Looks like we need to do a glue!
			glueResult = GlueEntry();
			glueSingleEntry(e, newEntry, key_norm, glueResult); //glue the two strings
			glueStorage.insertAfterFind(GlueEntry()); //clear entry
			//glueStorage.insertAtKey(key_norm, GlueEntry()); //clear entry
		}
	}

    nbGlueInserts++;
    // update stats and cleanup queue every 10M inserts
    if (nbGlueInserts % 1000000 == 0)
    {
        glueStorage.updateMemStats();
        glue();
    } 

	return true; 

	//if (key == key_of_interest) cout << "after\t" << glueStorage.dump(key,false ) << endl;
	if (glueDebug) cout << "after\t" << glueStorage.dump(key_norm,false ) << endl;

}


void GlueStorage::cleanup() {
	//GlueEntry e("bla", false, true, kmerSize);
	GlueEntry e;
	for (auto it = glueMap.begin(); it != glueMap.end(); ) {
        // later, TODO acquire a lock here so that it can be called safely in threads
		derefIt(it, e);
		if (e.isEmpty()) {
#ifdef SPARSEHASH
			glueMap.erase(it);
#else
			it = glueMap.erase(it);
#endif
		} else {
			it++;
		}
	}
#ifdef SPARSEHASH
	glueMap.resize(0); // effectively remove erased entries from memory
#endif
}


void Glue::glue()
{
	startTimer();
	glueStorage.cleanup(); 
	stopTimer();
}

void GlueCommander::output(string seq)
{
	Sequence s (Data::ASCII);
	s.getData().setRef ((char*)seq.c_str(), seq.size());
	out->insert(s);
}

void Glue::stopTimer() { 
	if (--timerReferenceCount == 0) {
		auto endTime = chrono::system_clock::now();
		totalTime += chrono::duration_cast<chrono::microseconds>(endTime - startTime);
	}
}

unsigned long GlueStorage::glueMapSize() {
    return glueMap.size();
}

unsigned long GlueStorage::glueMapSizeNonEmpty() {
    unsigned long res = 0;
	GlueEntry e;
	for (auto it = glueMap.begin(); it != glueMap.end(); ) {
		derefIt(it, e);
		if (!e.isEmpty()) {
            res++;
		}
		it++;
	}
    return res;
}



void GlueStorage::updateMemStats() {
	maxEntries = std::max(maxEntries, glueMap.size());
	totEntries += glueMap.size();

	size_t size = 0;
	for (auto it_gS = glueMap.begin(); it_gS != glueMap.end(); it_gS++)
	{
		size += sizeof(it_gS->first) + sizeof(it_gS->second);
		size += it_gS->first.size() + it_gS->second.size();
	}

	maxSize = std::max(maxSize, size);
	totSize += size;
	numDataPoints++;
}

void GlueStorage::printMemStats() {
	if (numDataPoints == 0) {
		cout << "GlueStorage: no data points to output memory stats.\n";
	} else {
		cout << "GlueStorage memory stats: max entries: " << add_commas(maxEntries) << " avg entries: " << add_commas(totEntries/numDataPoints) << " max size: " << add_commas(maxSize) 
			<< "b avg size: " << add_commas(totSize / numDataPoints) << "b\n";
	}
}



/*
 *
 * GlueCommander
 *
 */

void GlueCommander::cleanup()
{
    for (int i = 0; i < nb_glues; i++)
        glues[i]->glueStorage.cleanup();
}

void GlueCommander::updateMemStats()
{
    for (int i = 0; i < nb_glues; i++)
        glues[i]->glueStorage.updateMemStats();
}

void GlueCommander::dump()
{
    for (int i = 0; i < nb_glues; i++)
    {
        cout << "Printing memory stats for glue " << i << "\n";
        cout << glues[i]->glueStorage.dump();
    }
}

void GlueCommander::stop()
{
    cout << "Stopping GlueCommander" << endl;

    while (queues_size(true) != 0)
    {
        // wait for all queues and glues to be empty
        sleep(0.5);
    }

    for (int i = 0; i < nb_glues; i++)
        keep_glueing[i] = false;
    
    for (int i = 0; i < nb_glues; i++)
        glue_threads[i]->join(); 

    cleanup();    
}


void GlueCommander::printMemStats()
{
    for (int i = 0; i < nb_glues; i++)
    {
        cout << "Printing memory stats for glue " << i << "\n";
        glues[i]->glueStorage.printMemStats();
    }
}

int GlueCommander::which_queue(size_t minimizer)
{
    return minimizer % nb_glues;
}

void GlueCommander::spawn_threads()
{

    for (int i = 0; i < nb_glues; i++)
    {
        auto lambdaGlue = [this, i]() {

            std::pair<GlueEntry, string> glue_elt;
            cout << "Glue thread "<< i << " START" << endl;
            while (keep_glueing[i])
            {
                while (insert_aux_queues[i].try_dequeue(glue_elt))
                {
                    GlueEntry insRes;
                    glues[i]->insert_aux(glue_elt.first, glue_elt.second, insRes);
                    if (!insRes.isEmpty()) insert(insRes);
                }
            }
            cout << "Glue thread " << i << "END" << endl;
        };

        glue_threads.push_back(new std::thread (lambdaGlue));
    }

}

GlueCommander::GlueCommander(size_t _kmerSize, BankFasta *out, int nb_glues, Model *model) : nb_glues(nb_glues), model(model), out(out)
{

    for (int i = 0; i < nb_glues; i++)
    {
        glues.push_back(new Glue(_kmerSize, this));
        keep_glueing.push_back(true);
    }

    insert_aux_queues.resize(nb_glues);
    spawn_threads();
}


unsigned long GlueCommander::queues_size(bool silent)
{
    unsigned long sizes = 0;
    for (int i = 0; i < nb_glues; i++)
    {
        //if (!silent)
        {
            cout << "Size of insert_aux queue for glue " << i << " : " << insert_aux_queues[i].size_approx() << endl;
            cout << "Size of glueMap for glue " << i << " : " << glues[i]->glueStorage.glueMapSize() << endl;
        }
        //else
            sizes += insert_aux_queues[i].size_approx() + glues[i]->glueStorage.glueMapSizeNonEmpty() ;
    }
    return sizes;
}



void GlueCommander::insert(GlueEntry &e) {
	
	/* for testing:
	string str;
	while (cin >> str) {
		GlueEntry en(str, true, true, 2);
		RawEntry raw = en.getRaw();
		cout << tostring(en, "Z") << endl << raw << endl << tostring(GlueEntry(raw,2), "Z") << endl;
	}
	return; 	
	*/

	bool oldGlueDebug = glueDebug;
	
	//if ((e.getSeq().find(key_of_interest) != std::string::npos) || (e.getSeq().find(reversecomplement(key_of_interest)) != std::string::npos))  glueDebug = true; 

	if (glueDebug) 
		cout << "insert\t" << tostring(e,"") << endl;

	if (!e.getRmark() && !e.getLmark()) {
        commander_mutex.lock(); // hopefully not too much performance issue? else use a queue..
		output(e.getSeq());
        commander_mutex.unlock();
	} else if (e.getLmark() && !e.getRmark()) {
		insert_aux(e, e.getLkmer());
	} else if (!e.getLmark() && e.getRmark()) {
		insert_aux(e, e.getRkmer());
	} else { 
		//pick one that is not empty
		//this is not necessary for correctness, picking an arbitrary one will do
		//but picking a non-empty one may force glues to happen sooner rather later, preventing build-ups of long sequential chains
		/*if (!insert_aux(e, e.getRkmer(), insRes, true)) {
			insert_aux(e, e.getLkmer(), insRes);
		}*/

        // but for threaded glues, it's not so easy, to let's pick an arbitrary one
		insert_aux(e, e.getLkmer());
	}

	glueDebug = oldGlueDebug; 
}

bool GlueCommander::insert_aux(GlueEntry newEntry, string key) { 


    // maybe try to move that into Glue::insert_aux again but i'm not sure if it affects minimizer computation so i'm keeping it here to make sure 
	string key_norm = rcnorm(key);
	if (key != key_norm) {
		newEntry.revComp();
	}

    Model::Kmer kmer= model->codeSeed(key_norm.c_str(), Data::ASCII);
    size_t minimizer(model->getMinimizerValue(kmer.value()));
    int w = which_queue(minimizer); 
    insert_aux_queues[w].enqueue(make_pair(newEntry, key_norm));
}
