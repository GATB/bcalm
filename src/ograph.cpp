#include "ograph.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include "binSeq.h"
#include <thread>
#include <list> /* list */

#define BASE 256 // # of buckets to use

/*
 * constructs the compressed dBG from a list of arbitrary sequences
 */

using namespace std;

#ifdef SPARSE_HASH
	using namespace google;
#endif /* SPARSE_HASH */


struct comparator{bool operator()(const kmerIndice& a , const kmerIndice& b) { return a.kmmer < b.kmmer; }};


struct comparator2{bool operator()(const kmer2Indice& a , const kmer2Indice& b) { return min(a.indiceL,a.indiceR) < min(b.indiceL,b.indiceR); }};


struct equalo{
	bool operator()(const kmerIndice& a , const kmerIndice& b) { return a.kmmer == b.kmmer; }
};


bool compare (kmerIndice i,kmerIndice j) { return (i.kmmer<j.kmmer); }


void radix(vector<kmerIndice>& nums, uint64_t max) {
	vector<kmerIndice> bucket[BASE];
	uint64_t i;

	// iterate through each radix until n>max
	for (__uint128_t n=1; max >= n; n *= BASE) {
		// sort list of numbers into buckets
		for (i=0; i<nums.size(); i++){
			bucket[(nums[i].kmmer/n)%BASE].push_back(nums[i]);
		}

		// merge buckets back to list
		for (uint64_t k=i=0; i<BASE; bucket[i++].clear()){
			for (vector<kmerIndice>::iterator j = bucket[i].begin(); j != bucket[i].end(); nums[k++] = *(j++));
		}
	}
}


void radixkmer2indice(vector<kmer2Indice>& nums, uint64_t max) {
	vector<kmer2Indice> bucket[BASE];
	uint64_t i;

	// iterate through each radix until n>max
	for (__uint128_t n=1; max >= n; n *= BASE) {
		// sort list of numbers into buckets
		for (i=0; i<nums.size(); i++){
			bucket[(nums[i].indiceL/n)%BASE].push_back(nums[i]);
		}

		// merge buckets back to list
		for (uint64_t k=i=0; i<BASE; bucket[i++].clear()){
			for (vector<kmer2Indice>::iterator j = bucket[i].begin(); j != bucket[i].end(); nums[k++] = *(j++));
		}
	}
}


uint nt2num(char c){
	// 'a' -> 0, 'c' -> 1; 'g' -> 2, 't' -> 3
	// inspired by G. Rizk
	char d = (c >> 1) & 3;
	if (d > 1)
		d ^= 1;
	return d;
}


char num2nt(int num){
	if (num == 0)
		return 'a';
	else if (num == 1)
		return 'c';
	else if (num == 2)
		return 'g';
	else if (num == 3)
		return 't';
	assert(0);
	return '*';
}


uint nt2num(char c, int pos){
	int offset = pos % 4;
	return (nt2num(c) + offset ) % 4;
}


char num2nt(int num, int pos){
	int offset = pos % 4;
	return num2nt((num - offset + 4) % 4);
}


int minbutbiggerthan(int m1, int m2, const string &namebucket){
	int h(stoi(namebucket));
	if(m1<m2){
		if(m1>h)
			return m1;
		if(m2>h)
			return m2;
	}
	else{
		if(m2>h)
			return m2;
		if(m1>h)
			return m1;
	}
	return -1;
}


string reversecompletment(const string& str){
	string res(str);
	for(int i(str.size()-1), j(0); i > -1; --i, ++j){
		unsigned char c = str[i];
		c ^= 4;
		if ((c&3) != 3){c ^= 17;}
		res[j] = c;
		// res[i]=3-nt2num(str[j]);
	}
	return res;
}


string reverseinplace(string& str){
	uint i(str.size()-1),j(0);
	for(; j<str.size()/2; --i, ++j){
		// unsigned char c = str[i];
		// unsigned char d = str[j];
		// c ^= 4;
		// d ^= 4;
		// if ((c&3) != 3){c ^= 17;}
		// if ((d&3) != 3){d ^= 17;}
		// str[j]=c;
		// str[i]=d;

		str[i] ^= 4;
		str[j] ^= 4;
		if ((str[i]&3) != 3){str[i]^= 17;}
		if ((str[j]&3) != 3){str[j]^= 17;}
		swap(str[i],str[j]);
	}
	if(str.size()%2==1){
		str[j] ^= 4;
		if ((str[j]&3) != 3){str[j]^= 17;}
	}
	return str;
}


void reverseinplace2(string& str){
	uint i(str.size()-1),j(0);
	for(; j<str.size()/2; --i, ++j){
		// unsigned char c = str[i];
		// unsigned char d = str[j];
		// c ^= 4;
		// d ^= 4;
		// if ((c&3) != 3){c ^= 17;}
		// if ((d&3) != 3){d ^= 17;}
		// str[j]=c;
		// str[i]=d;

		str[i] ^= 4;
		str[j] ^= 4;
		if ((str[i]&3) != 3){str[i]^= 17;}
		if ((str[j]&3) != 3){str[j]^= 17;}
		swap(str[i],str[j]);
	}
	if(str.size()%2==1){
		str[j] ^= 4;
		if ((str[j]&3) != 3){str[j]^= 17;}
	}
}


string getRepresent(const string& str){
	return(min(str,reversecompletment(str)));
}


string readn(ifstream *file,uint64_t n){
	string contents(n,'0');
	file->read(&contents[0],n);
	return(contents);
}


bool adjacent(const string& node1,const  string& node2,int k){
	return(node1.substr(node1.size()-k+1,k-1)==node2.substr(0,k-1));
}


uint chartoint(char c){
	char d = (c >> 1) & 3;
	if (d > 1)
		d ^= 1;
	return d;
}


char complement(char c){
	switch(c){
		case 'a':
		return 't';
		case 'c':
		return 'g';
		case 'g':
		return 'c';
		case 't':
		return 'a';
		default:
		return 0;
	}
}


uint64_t stringtoint(const string& str){
	uint64_t res(0);
	for(uint32_t i(0);i<str.size();i++){
		res<<=2;
		res+=chartoint(str[i]);
	}
	return res;
}


uint64_t stringtointc(const string& str){
	uint64_t res(0);
	for(int32_t i(str.size()-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[i]);
	}
	return res;
}


__uint128_t stringtoint128(const string& str){
	__uint128_t res(0);
	for(uint32_t i(0);i<str.size();i++){
		res<<=2;
		res+=chartoint(str[i]);
	}
	return res;
}


__uint128_t graph3::beg2int128(const string& str){
	__uint128_t resBeg(0);
	for(uint i(0);i<k;++i){
		resBeg<<=2;
		resBeg+=chartoint(str[i]);
	}
	return resBeg;
}


__uint128_t graph3::beg2int128rc(const string& str){
	__uint128_t res(0);
	for(int32_t i(k-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[i]);
	}
	return res;
}


__uint128_t graph3::end2int128rc(const string& str){
	__uint128_t res(0);
	for(int32_t i(k-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[str.size()-k+i]);
	}
	return res;
}


__uint128_t graph3::end2int128(const string& str){
	__uint128_t resEnd(0);
	for(uint i(0);i<k;++i){
		resEnd<<=2;
		resEnd+=chartoint(str[str.size()-k+i]);
	}
	return resEnd;
}


__uint128_t stringtointc128(const string& str){
	__uint128_t res(0);
	for(int32_t i(str.size()-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[i]);
	}
	return res;
}


uint64_t string2intmin(const string& str){
	return min(stringtoint(str),stringtointc(str));
}


bool accordtomin(int min, int left_or_right_min){
	if(min == -1){
		return true;
	}

	if(left_or_right_min==min)
		return true;

	return false;

}


string compaction2(const string& seq1,const string& seq2, int k){
	size_t s1(seq1.size()),s2(seq2.size());
	if(s1==0){return seq2;}
	if(s2==0){return seq1;}

	string rc2(reversecompletment(seq2));
	string rc1(reversecompletment(seq1));


	if(seq1.substr(0,k)==seq2.substr(s2-k,k)){
		return seq2+seq1.substr(k);
	}else{
		if(rc2.substr(s2-k,k)==seq1.substr(0,k)){
			return rc2+seq1.substr(k);
		}
	}

	if(seq2.substr(0,k)==seq1.substr(s1-k,k)){
		return seq1+seq2.substr(k);
	}else{
		if(rc1.substr(s1-k,k)==seq2.substr(0,k)){
			return rc1+seq2.substr(k);
		}
	}

	if(rc1.substr(0,k)==seq2.substr(s2-k,k)){
			return seq2+rc1.substr(k);
	}else{
		if(rc2.substr(s2-k,k)==rc1.substr(0,k)){
			return rc2+rc1.substr(k);
		}
	}

	if(rc2.substr(0,k)==seq1.substr(s1-k,k)){
		return seq1+rc2.substr(k);
	}else{
		if(rc1.substr(s1-k,k)==rc2.substr(0,k)){
			return rc1+rc2.substr(k);
		}
	}
	return seq1;
}


string compaction(const string& seq1,const string& seq2, int k){
	int s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){
		return seq1;
	}
	string beg1(seq1.substr(0,k));
	if(beg1==seq2.substr(s2-k,k)){
		return seq2+seq1.substr(k);
	}

	string end1(seq1.substr(s1-k,k));
	if(seq2.substr(0,k)==end1){
		return seq1+seq2.substr(k);
	}

	string rc2(reversecompletment(seq2));
	if(rc2.substr(s2-k,k)==beg1){
		return rc2+seq1.substr(k);
	}

	if(rc2.substr(0,k)==end1){
		return seq1+rc2.substr(k);
	}

	return seq1;
}


string compactionBeg(const string& seq1,const string& seq2, int k){
	int s2(seq2.size());
	string rc2(reversecompletment(seq2));
	//~ string rc1(reversecompletment(seq1));
	string beg(seq1.substr(0,k));
	string begRc(reversecompletment(beg));


	if(beg==seq2.substr(s2-k,k)){
		return seq2+seq1.substr(k);
	}else{
		if(beg==rc2.substr(s2-k,k)){
			return rc2+seq1.substr(k);
		}
	}

	if(begRc==seq2.substr(0,k)){
		return rc2.substr(0,s2-k)+seq1;
	}else{
		if(begRc==rc2.substr(0,k)){
			return seq2.substr(0,s2-k)+seq1;
		}
	}
	return seq1;
}


string compactionEnd(const string& seq1,const string& seq2, int k){
	int s1(seq1.size()),s2(seq2.size());
	string rc2(reversecompletment(seq2));
	string end(seq1.substr(s1-k,k));
	string endRc(reversecompletment(end));


	if(end==seq2.substr(0,k)){
		return seq1+seq2.substr(k);
	}else{
		if(end==rc2.substr(0,k)){
			return seq1+rc2.substr(k);
		}
	}

	if(endRc==seq2.substr(s2-k,k)){
		return seq1+rc2.substr(k);
	}else{
		if(endRc==rc2.substr(s2-k,k)){
			return seq1+seq2.substr(k);
		}
	}
	return seq1;
}


bool isNumber(char c){return (c<64);}


__uint128_t graph3::rcb(__uint128_t min){
	__uint128_t resrcb(0);
	__uint128_t offsetrcb(1);
	// offsetrcb<<=(2*k-2);
	for(uint i(0); i<k;++i){
		resrcb+=(3-(min%4))<<(2*(k-1-i));
		min>>=2;
		// offsetrcb>>=2;
	}
	return resrcb;
}


void graph3::compaction2( uint32_t iL,  uint32_t iR){
	uint s1(unitigs[iL].size()),s2(unitigs[iR].size());
	bool b1(isNumber(unitigs[iL][0])),b2(isNumber(unitigs[iR][0]));
	if(b1 and b2){return compaction(stoi(unitigs[iL]),stoi(unitigs[iR]));}
	if(b1){return compaction(stoi(unitigs[iL]),iR);}
	if(b2){return compaction(iL,stoi(unitigs[iR]));}
	__uint128_t beg2(beg2int128(unitigs[iR]));
	__uint128_t end1(end2int128(unitigs[iL]));
	if(end1==beg2){
		++u1;
		unitigs[iL]+=(unitigs[iR].substr(k));
		unitigs[iR]=to_string(iL);
		return;
	}
	__uint128_t begrc2(end2int128rc(unitigs[iR]));
	if(end1==begrc2){
		++u2;
		unitigs[iL]+=(reverseinplace(unitigs[iR]).substr(k));
		unitigs[iR]=to_string(iL);
		return;
	}
	__uint128_t beg1(beg2int128(unitigs[iL]));
	__uint128_t end2(rcb(begrc2));
	if(beg1==end2){
		++u3;
		unitigs[iR]+=(unitigs[iL].substr(k));
		unitigs[iL]=to_string(iR);
		return;
	}
	__uint128_t endrc2(rcb(beg2));
	if(beg1==endrc2){
		++u4;
		reverseinplace2(unitigs[iR]);
		unitigs[iR]+=(unitigs[iL].substr(k));
		unitigs[iL]=to_string(iR);
		return;
	}
}


void graph3::compaction( uint32_t iL,  uint32_t iR){
	uint s1(unitigs[iL].size()),s2(unitigs[iR].size());
	bool b1(isNumber(unitigs[iL][0])),b2(isNumber(unitigs[iR][0]));
	if(b1 and b2){return compaction(stoi(unitigs[iL]),stoi(unitigs[iR]));}
	if(b1){return compaction(stoi(unitigs[iL]),iR);}
	if(b2){return compaction(iL,stoi(unitigs[iR]));}

	__uint128_t beg1(beg2int128(unitigs[iL]));
	__uint128_t end2(end2int128(unitigs[iR]));
	if(beg1==end2){
		unitigs[iR]+=(unitigs[iL].substr(k));
		unitigs[iL]=to_string(iR);
		return;
	}

	__uint128_t endrc2(beg2int128rc(unitigs[iR]));
	if(beg1==endrc2){
		reverseinplace2(unitigs[iR]);
		unitigs[iR]+=(unitigs[iL].substr(k));
		unitigs[iL]=to_string(iR);
		return;
	}

	__uint128_t beg2(rcb(endrc2));
	__uint128_t end1(end2int128(unitigs[iL]));
	if(end1==beg2){
		unitigs[iL]+=(unitigs[iR].substr(k));
		unitigs[iR]=to_string(iL);
		return;
	}

	__uint128_t begrc2(rcb(end2));
	if(end1==begrc2){
		unitigs[iL]+=(reverseinplace(unitigs[iR]).substr(k));
		unitigs[iR]=to_string(iL);
		return;
	}
}


void makeUnique(vector<kmerIndice>& V){
	uint j(0);
	bool take(true);
	kmerIndice previous(V[0]);
	for(uint i(1);i<V.size();++i){
		if(previous.kmmer==V[i].kmmer){
			take=false;
		}else{
			if(take){
				V[j++]=previous;
			}else{
				take=true;
			}
			previous=V[i];
		}
	}
	V.resize(j);

}


void graph3::debruijn2(){
	sort(left.begin(),left.end(),comparator());
	sort(right.begin(),right.end(),comparator());
	makeUnique(left);
	makeUnique(right);
	uint iL(0),iR(0),sL(left.size()),sR(right.size());
	while(sL!=iL and sR!=iR){
		if(left[iL].kmmer==right[iR].kmmer){
			++iL;++iR;
			compaction(left[iL].indice,right[iR].indice);
		}else{
			if(left[iL].kmmer<right[iR].kmmer){
				++iL;
			}else{
				++iR;
			}
		}
	}
}


void graph3::debruijn(){
	sort(left.begin(),left.end(),comparator());
	sort(right.begin(),right.end(),comparator());
	uint iL(0),iR(0),sL(left.size()),sR(right.size());
	kmerIndice kL,kR;
	while(sL!=iL and sR!=iR){
		kL=left[iL];
		kR=right[iR];
		if(kL.kmmer==kR.kmmer){
			bool go(true);
			++iL;++iR;
			if(sL!=iL){
				if(left[iL].kmmer==kL.kmmer){
					++iL;
					go=false;
					while(sL!=iL){
						if(left[iL].kmmer!=kL.kmmer){break;}else{++iL;}
					}
				}
			}
			if(iR!=sR){
				if(right[iR].kmmer==kL.kmmer){
					++iR;
					go=false;
					while(sR!=iR){
						if(right[iR].kmmer!=kR.kmmer){break;}else{++iR;}
					}
				}
			}
			if(go){
				compaction(kL.indice,kR.indice);
			}
		}else{
			if(kL.kmmer<kR.kmmer){
				++iL;
				while(sL!=iL){
					if(left[iL].kmmer!=kL.kmmer){break;}else{++iL;}
				}
			}else{
				++iR;
				while(sR!=iR){
					if(right[iR].kmmer!=kR.kmmer){break;}else{++iR;}
				}
			}
		}
	}
}


bool graph3::output(uint i){
return !isNumber(unitigs[i][0]);
}


bool graph3::clear(){delete [] unitigs;return true;}


uint32_t graph3::size(){return indiceUnitigs;};


void graph3::addtuple(tuple<string,uint32_t,uint32_t>& tuple){
	unitigs[indiceUnitigs]=move(get<0>(tuple));
	if(minimizer==(get<1>(tuple))){
		__uint128_t kmer1(beg2int128(unitigs[indiceUnitigs]));
		// kmer2=(beg2int128rc(unitigs[indiceUnitigs]));
		__uint128_t kmer2(rcb(kmer1));
		// ki.indice=indiceUnitigs;
		if(kmer1<kmer2){
			left.push_back(kmerIndice{indiceUnitigs,kmer1});
		}else{
			right.push_back(kmerIndice{indiceUnitigs,kmer2});
		}
	}
	if(minimizer==get<2>(tuple)){
		__uint128_t kmer1(end2int128(unitigs[indiceUnitigs]));
		// kmer2=(end2int128rc(unitigs[indiceUnitigs]));
		__uint128_t kmer2(rcb(kmer1));
		// ki.indice=indiceUnitigs;
		if(kmer1<kmer2){
			// ki.kmmer=kmer1;
			right.push_back(kmerIndice{indiceUnitigs,kmer1});
		}else{
			// ki.kmmer=kmer2;
			left.push_back(kmerIndice{indiceUnitigs,kmer2});
		}
	}
	++indiceUnitigs;
}


void TS(vector<kmerIndice>* V){
	sort(V->begin(),V->end(),comparator());
}

