#include "ograph.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include "binSeq.h"


/*
 * constructs the compressed dBG from a list of arbitrary sequences
 */


using namespace std;


struct comparator{bool operator()(const kmerIndice& a , const kmerIndice& b) { return a.kmmer < b.kmmer; }};


string reverseinplace(string& str){
	uint i(str.size()-1),j(0);
	for(; j<str.size()/2; --i, ++j){
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


uint chartoint(char c){
	char d = (c >> 1) & 3;
	if (d > 1)
		d ^= 1;
	return d;
}


bool isNumber(char c){return (c<64);}


kmer graph4::beg2int128(const string& str){
	kmer resBeg(0);
	for(uint i(0);i<k;++i){
		resBeg<<=2;
		resBeg+=chartoint(str[i]);
	}
	return resBeg;
}


kmer graph4::beg2int128rc(const string& str){
	kmer res(0);
	for(int i(k-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[i]);
	}
	return res;
}


kmer graph4::end2int128rc(const string& str){
	kmer res(0);
	for(int i(k-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[str.size()-k+i]);
	}
	return res;
}


kmer graph4::end2int128(const string& str){
	kmer resEnd(0);
	for(uint i(0);i<k;++i){
		resEnd<<=2;
		resEnd+=chartoint(str[str.size()-k+i]);
	}
	return resEnd;
}


kmer graph4::rcb(kmer min){
	kmer resrcb(0);
	kmer offsetrcb(1);
	for(uint i(0); i<k;++i){
		resrcb+=(3-(min%4))<<(2*(k-1-i));
		min>>=2;
	}
	return resrcb;
}


void graph4::compaction( uint iL,  uint iR){
	// cout<<iL<<" "<<iR<<endl;
	uint s1(unitigs[iL].size()),s2(unitigs[iR].size());
	bool b1(unitigs[iL].isInt),b2(unitigs[iR].isInt);
	if(b1 and b2){return compaction((unitigs[iL].getNumber()),(unitigs[iR]).getNumber());}
	if(b1){return compaction(unitigs[iL].getNumber(),iR);}
	if(b2){return compaction(iL,unitigs[iR].getNumber());}

	kmer end1(unitigs[iL].getEndInt(k));
	kmer beg2(unitigs[iR].getBeginInt(k));
	if(end1==beg2){
		unitigs[iL].add(unitigs[iR].sub(k));
		unitigs[iR]=binSeq(iL);
		return;
	}

	binSeq rc2(unitigs[iR].getReverse());
	kmer begrc2(unitigs[iR].getEndRcInt(k));
	//~ __uint128_t begrc2(rc2.getBeginInt(k));
	if(end1==begrc2){
		unitigs[iL].add(rc2.sub(k));
		unitigs[iR]=binSeq(iL);
		return;
	}

	kmer beg1(unitigs[iL].getBeginInt(k));
	kmer end2(unitigs[iR].getEndInt(k));
	if(end2==beg1){
		unitigs[iR].add(unitigs[iL].sub(k));
		unitigs[iL]=binSeq(iR);
		return;
	}

	//~ __uint128_t endrc2(rc2.getEndInt(k));
	kmer endrc2(unitigs[iR].getBeginRcInt(k));
	if(endrc2==beg1){
		rc2.add(unitigs[iL].sub(k));
		unitigs[iR]=rc2;
		unitigs[iL]=binSeq(iR);
		return;
	}
// cout<<"wut"<<endl;
}


void graph4::debruijn(){
	sort(left.begin(),left.end(),comparator());
	sort(right.begin(),right.end(),comparator());
	uint iL(0),iR(0);
	kmerIndice kL,kR;
	while(iL!=left.size() and iR!=right.size()){
		kL=left[iL];
		kR=right[iR];
		if(kL.kmmer==kR.kmmer){
			bool go(true);
			++iL;++iR;
			if(left[iL].kmmer==kL.kmmer){
				go=false;
				while(left[++iL].kmmer==kL.kmmer){}
			}
			if(right[iR].kmmer==kL.kmmer){
				go=false;
				while(right[++iR].kmmer==kR.kmmer){}
			}
			if(go){compaction(kL.indice,kR.indice);}
		}else{
			if(kL.kmmer<kR.kmmer){
				while(++iL!=left.size() and left[iL].kmmer<kR.kmmer){}
			}else{
				while(++iR!=right.size() and right[iR].kmmer<kL.kmmer){}
			}
		}
	}
}


bool graph4::output(uint i){return (!unitigs[i].isInt);}


void graph4::clear(){delete [] unitigs;}


uint graph4::size(){return indiceUnitigs;};


void graph4::addtuple(tuple<binSeq,uint,uint>& tuple){
	unitigs[indiceUnitigs]=(get<0>(tuple));
	if(minimizer==(get<1>(tuple))){
		kmer kmer1(unitigs[indiceUnitigs].getBeginInt(k));
		//~ beg.reverse();
		kmer kmer2(unitigs[indiceUnitigs].getBeginRcInt(k));
		// kmer kmer1(beg2int128(unitigs[indiceUnitigs]));
		// kmer kmer2(rcb(kmer1));
		if(kmer1<kmer2){
			left.push_back(kmerIndice{indiceUnitigs,kmer1});
		}else{
			right.push_back(kmerIndice{indiceUnitigs,kmer2});
		}
	}
	if(minimizer==get<2>(tuple)){
		kmer kmer1(unitigs[indiceUnitigs].getEndInt(k));
		//~ end.reverse();
		kmer kmer2(unitigs[indiceUnitigs].getEndRcInt(k));
		// kmer kmer1(end2int128(unitigs[indiceUnitigs]));
		// kmer kmer2(rcb(kmer1));
		if(kmer1<kmer2){
			right.push_back(kmerIndice{indiceUnitigs,kmer1});
		}else{
			left.push_back(kmerIndice{indiceUnitigs,kmer2});
		}
	}
	++indiceUnitigs;
}
