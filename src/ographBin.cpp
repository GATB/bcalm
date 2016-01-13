#include "ographBin.h"
#include "ograph.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include "binSeq.h"


/*
 * constructs the compressed dBG from a list of arbitrary sequences
 */


using namespace std;


struct comparatorOver{bool operator()(const overlapIndice& a , const overlapIndice& b) { return a.kmmer < b.kmmer; }};


kmer graph4::rcb(kmer min){
	kmer resrcb(0);
	kmer offsetrcb(1);
	for(uint i(0); i<k;++i){
		resrcb+=(3-(min%4))<<(2*(k-1-i));
		min>>=2;
	}
	return resrcb;
}


void graph4::compaction(uint iL,  uint iR){
	// cout<<iL<<" "<<iR<<endl;
	// cout<<"gocomp"<<endl;
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
	// cout<<"debruijn go"<<endl;
	sort(left.begin(),left.end(),comparatorOver());
	sort(right.begin(),right.end(),comparatorOver());
	uint iL(0),iR(0);
	overlapIndice kL,kR;
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
	// cout<<"debruijn end"<<endl;
}


bool graph4::output(uint i){return (!unitigs[i].isInt);}


void graph4::clear(){delete [] unitigs;}


uint graph4::size(){return indiceUnitigs;};


void graph4::addtuple(tuple<binSeq,uint,uint>& tuple){
	unitigs[indiceUnitigs]=(get<0>(tuple));
	if(minimizer==(get<1>(tuple))){
		overlap kmer1(unitigs[indiceUnitigs].getBeginIntOver(k));
		overlap kmer2(unitigs[indiceUnitigs].getBeginRcIntOver(k));
		if(kmer1<kmer2){
			left.push_back(overlapIndice{indiceUnitigs,kmer1});
		}else{
			right.push_back(overlapIndice{indiceUnitigs,kmer2});
		}
	}
	if(minimizer==get<2>(tuple)){
		overlap kmer1(unitigs[indiceUnitigs].getEndIntOver(k));
		overlap kmer2(unitigs[indiceUnitigs].getEndRcIntOver(k));
		if(kmer1<kmer2){
			right.push_back(overlapIndice{indiceUnitigs,kmer1});
		}else{
			left.push_back(overlapIndice{indiceUnitigs,kmer2});
		}
	}
	++indiceUnitigs;
}
