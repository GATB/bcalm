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
	int n = str.size();
	for(int i(n-1), j(0); i > -1; i--, j++){
		unsigned char c = str[i];
		unsigned char d = (c >> 4)&7;
		c ^= 4;
		if ((c&3) != 3)
			c ^= 17;
		res[j] = c;
		}
	return res;
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
	__uint128_t resBeg=(0);
	for(uint i(0);i<k;i++){
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
	__uint128_t resEnd=(0);
	for(uint i(0);i<k;i++){
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
	//~ cout<<"beg"<<endl;
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
		//~ cout<<"beg"<<endl;
	//~ cout<<seq1<<" "<<seq2<<endl;


	return seq1;
}


string compactionEnd(const string& seq1,const string& seq2, int k){
	int s1(seq1.size()),s2(seq2.size());
	string rc2(reversecompletment(seq2));
	//~ string rc1(reversecompletment(seq1));
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
	//~ cout<<"end"<<endl;
	//~ cout<<seq1<<" "<<seq2<<endl;

	return seq1;
}


bool isNumber(char c){
	return (c<64);
}


void graph3::compaction2(uint32_t iL, uint32_t iR){
	string seq1(unitigs[iL]),seq2(unitigs[iR]);
	size_t s1(seq1.size()),s2(seq2.size());

	bool b1(isNumber(seq1[0])),b2(isNumber(seq2[0]));
	if(b1 && b2){
		compaction(stoi(seq1),stoi(seq2));
		return;
	}
	if(b1){
		compaction(stoi(seq1),iR);
		return;
	}
	if(b2){
		compaction(iL,stoi(seq2));
		return;
	}

	string rc2(reversecompletment(seq2));
	string rc1(reversecompletment(seq1));

	if(seq1.substr(0,k)==seq2.substr(s2-k,k)){
		// cout<<1<<endl;
		unitigs[iL]=seq2+seq1.substr(k);
		unitigs[iR]=to_string(iL);
		return;
	}else{
		if(rc2.substr(s2-k,k)==seq1.substr(0,k)){
				// cout<<2<<endl;
		unitigs[iL]=rc2+seq1.substr(k);
		unitigs[iR]=to_string(iL);
		return;
		}
	}
	if(seq2.substr(0,k)==seq1.substr(s1-k,k)){
			// cout<<3<<endl;
		unitigs[iL]=seq1+seq2.substr(k);
		unitigs[iR]=to_string(iL);
		return;
	}else{
		if(rc1.substr(s1-k,k)==seq2.substr(0,k)){
				// cout<<4<<endl;
			unitigs[iL]=rc1+seq2.substr(k);
			unitigs[iR]=to_string(iL);
			return;
		}
	}
	if(rc1.substr(0,k)==seq2.substr(s2-k,k)){
			// cout<<5<<endl;
			unitigs[iL]=seq2+rc1.substr(k);
			unitigs[iR]=to_string(iL);
			return;
	}else{
		if(rc2.substr(s2-k,k)==rc1.substr(0,k)){
				// cout<<6<<endl;
			unitigs[iL]=rc2+rc1.substr(k);
			unitigs[iR]=to_string(iL);
			return;
		}
	}
	if(rc2.substr(0,k)==seq1.substr(s1-k,k)){
			// cout<<7<<endl;
		unitigs[iL]=seq1+rc2.substr(k);
		unitigs[iR]=to_string(iL);
		return;
	}else{
		if(rc1.substr(s1-k,k)==rc2.substr(0,k)){
				// cout<<8<<endl;
			unitigs[iL]=rc1+rc2.substr(k);
			unitigs[iR]=to_string(iL);
			return;
		}
	}
	cout<<"WUT"<<endl;
	cout<<"+"<<seq1<<"+ +"<<seq2<<"+"<<endl;
	cin.get();
}


__uint128_t graph3::rcb(__uint128_t min){
	__uint128_t resrcb(0);
	__uint128_t offsetrcb(1);
	offsetrcb<<=(2*k-2);
	for(uint i(0); i<k;++i){
		resrcb+=(3-(min%4))*offsetrcb;
		min>>=2;
		offsetrcb>>=2;
	}
	return resrcb;
}


void graph3::compaction( uint32_t iL,  uint32_t iR){
	uint s1(unitigs[iL].size()),s2(unitigs[iR].size());
	bool b1(isNumber(unitigs[iL][0])),b2(isNumber(unitigs[iR][0]));
	if(b1 and b2){return compaction(stoi(unitigs[iL]),stoi(unitigs[iR]));}
	if(b1){return compaction(stoi(unitigs[iL]),iR);}
	if(b2){return compaction(iL,stoi(unitigs[iR]));}
	__uint128_t beg2(beg2int128(unitigs[iR]));
	__uint128_t end1(end2int128(unitigs[iL]));
	// string end1(seq1.substr(s1-k,k));
	// string beg2(seq2.substr(0,k));
	if(end1==beg2){
	// if(true){
		unitigs[iL]+=(unitigs[iR].substr(k));
		unitigs[iR]=to_string(iL);
		return;
	}
	__uint128_t begrc2(end2int128rc(unitigs[iR]));
	// string rc2(reversecompletment(seq2));
	// string begrc2(rc2.substr(0,k));
	if(end1==begrc2){
		unitigs[iL]+=(reversecompletment(unitigs[iR]).substr(k));
		unitigs[iR]=to_string(iL);
		return;
	}
	__uint128_t beg1(beg2int128(unitigs[iL]));
	__uint128_t end2(rcb(begrc2));
	// string beg1(seq1.substr(0,k));
	// string end2(seq2.substr(s2-k,k));
	if(beg1==end2){
		unitigs[iR]+=(unitigs[iL].substr(k));
		unitigs[iL]=to_string(iR);
		return;
	}
	// endrc2=(beg2int128rc(seq2));
	__uint128_t endrc2(rcb(beg2));
	// string endrc2(rc2.substr(s2-k,k));
	if(beg1==endrc2){
		unitigs[iL]=reversecompletment(unitigs[iR])+(unitigs[iL].substr(k));
		unitigs[iR]=to_string(iL);
		return;
	}
	cout<<"wut"<<endl;
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
			}else{break;}
			if(iR!=sR){
				if(right[iR].kmmer==kL.kmmer){
					++iR;
					go=false;
					while(sR!=iR){
						if(right[iR].kmmer!=kR.kmmer){break;}else{++iR;}
					}
				}
			}else{break;}
			if(go){
				// kmer2Indice k2i;
				// k2i.indiceL=kL.indice;
				// k2i.indiceR=kR.indice;
				// compactions.push_back(k2i);
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
					if(right[iR++].kmmer!=kR.kmmer){break;}else{++iR;}
				}
			}
		}
	}
}


bool graph3::output(uint i){return !isNumber(unitigs[i][0]);}


bool graph3::clear(){delete [] unitigs;}


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


void graph4::debruijn(){
	sort(left.begin(),left.end(),comparator());
	sort(right.begin(),right.end(),comparator());
	isNumber.assign(unitigs.size(),false);
	kmerIndice kL,kR;
	while(left.size()!=0 and right.size()!=0){
		kL=left.back();
		kR=right.back();
		if(kL.kmmer==kR.kmmer){
			bool go(true);
			left.pop_back();
			right.pop_back();
			if(left.size()!=0){
				if(left.back().kmmer==kL.kmmer){
					left.pop_back();
					go=false;
					bool again;
					if(left.size()!=0){
						again=left.back().kmmer==kL.kmmer;
					}else{
						again=false;
					}
					while(again){
						left.pop_back();
						if(left.size()!=0){
							again=left.back().kmmer==kL.kmmer;
						}else{
							again=false;
						}
					}
				}
			}
			if(right.size()!=0){
				if(right.back().kmmer==kL.kmmer){
					right.pop_back();
					bool again;
					while(again){
						right.pop_back();
						if(right.size()!=0){
							again=right.back().kmmer==kR.kmmer;
						}else{
							again=false;
						}
					}
				}
			}
			if(go){
				//~ kmer2Indice k2i;
				//~ k2i.indiceL=min(kL.indice,kR.indice);
				//~ if(k2i.indiceR>maximum){maximum=k2i.indiceR;}
				//~ k2i.indiceR=max(kL.indice,kR.indice);
				//~ compactions.push_back(k2i);
				compaction(kL.indice,kR.indice);
			}
		}else{
			if(kL.kmmer>kR.kmmer){
				left.pop_back();
				bool again;
				if(left.size()!=0){
					again=left.back().kmmer==kL.kmmer;
				}else{
					again=false;
				}
				while(again){
					left.pop_back();
					if(left.size()!=0){
						again=left.back().kmmer==kL.kmmer;
					}else{
						again=false;
					}
				}
			}else{
				right.pop_back();
				bool again;
				if(right.size()!=0){
					again=right.back().kmmer==kL.kmmer;
				}else{
					again=false;
				}
				while(again){
					right.pop_back();
					if(right.size()!=0){
						again=right.back().kmmer==kR.kmmer;
					}else{
						again=false;
					}
				}
			}
		}
	}
}


void graph4::compress(){}


void graph4::addvertex(const string& unitigstr){
	//~ isNumber.push_back(false);
	binSeq unitig(unitigstr);
	unitigs.push_back(unitig);
	uint32_t i(unitigs.size()-1);

	if(leftmins[i]){
		//~ binSeq beg(unitig.getBegin(k));
		__uint128_t leftKmer1(unitig.getBeginInt(k));
		//~ beg.reverse();
		__uint128_t leftKmer2(unitig.getBeginRcInt(k));
		kmerIndice ki;
		ki.indice=i;

		if(leftKmer1<leftKmer2){
			ki.kmmer=leftKmer1;
			left.push_back(ki);
		}else{
			ki.kmmer=leftKmer2;
			right.push_back(ki);
		}
	}

	if(rightmins[i]){
		//~ binSeq end(unitig.getEnd(k));
		__uint128_t rightKmer1(unitig.getEndInt(k));
		//~ end.reverse();
		__uint128_t rightKmer2(unitig.getEndRcInt(k));
		kmerIndice ki;
		ki.indice=i;
		if(rightKmer1<rightKmer2){
			ki.kmmer=rightKmer1;
			right.push_back(ki);
		}else{
			ki.kmmer=rightKmer2;
			left.push_back(ki);
		}
	}
}


void graph4::addleftmin(unsigned int min){
	if(min==minimizer){
		leftmins.push_back(true);
	}else{
		leftmins.push_back(false);
	}
}


void graph4::addrightmin(unsigned int min){
	if(min==minimizer){
		rightmins.push_back(true);
	}else{
		rightmins.push_back(false);
	}
}


void print_bytes(void *p){
	size_t i;
	printf("(");
	for (i = 0; i < 8; ++i)
		printf("%02X", ((unsigned char*)p)[i]);
	printf(")");
}


void graph4::compaction(uint32_t iL, uint32_t iR){
	bool b1(isNumber[iL]),b2(isNumber[iR]);
	if(b1){
		if(b2){
			compaction(unitigs[iL].getNumber(),unitigs[iR].getNumber());
			return;
		}else{
			compaction(unitigs[iL].getNumber(),iR);
			return;
		}
	}else{
		if(b2){
			compaction(iL,unitigs[iR].getNumber());
			return;
		}
	}
		__uint128_t end1(unitigs[iL].getEndInt(k));
		__uint128_t beg2(unitigs[iR].getBeginInt(k));
		if(end1==beg2){
			unitigs[iL].add(unitigs[iR].sub(k));
			unitigs[iR]=binSeq(iL);
			isNumber[iR]=true;
			return;
		}


		binSeq rc2(unitigs[iR].getReverse());
		__uint128_t begrc2(unitigs[iR].getEndRcInt(k));
		//~ __uint128_t begrc2(rc2.getBeginInt(k));
		if(end1==begrc2){
			unitigs[iL].add(rc2.sub(k));
			unitigs[iR]=binSeq(iL);
			isNumber[iR]=true;
			return;
		}

		__uint128_t beg1(unitigs[iL].getBeginInt(k));
		__uint128_t end2(unitigs[iR].getEndInt(k));
		if(end2==beg1){
			unitigs[iR].add(unitigs[iL].sub(k));
			unitigs[iL]=binSeq(iR);
			isNumber[iL]=true;
			return;
		}

		//~ __uint128_t endrc2(rc2.getEndInt(k));
		__uint128_t endrc2(unitigs[iR].getBeginRcInt(k));
		if(endrc2==beg1){
			rc2.add(unitigs[iL].sub(k));
			unitigs[iR]=rc2;
			unitigs[iL]=binSeq(iR);
			isNumber[iL]=true;
			return;
		}
	cout<<"WUT"<<endl;
	cout<<"+"<<unitigs[iL].str()<<"+ +"<<unitigs[iR].str()<<"+"<<endl<<endl<<endl;
}
