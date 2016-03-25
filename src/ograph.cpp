#include "ograph.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <thread>


/*
 * constructs the compressed dBG from a list of arbitrary sequences
 */


using namespace std;


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


kmer graph3::beg2int128(const string& str){
	kmer resBeg(0);
	for(uint i(0);i<k;++i){
		resBeg<<=2;
		resBeg+=chartoint(str[i]);
	}
	return resBeg;
}


kmer graph3::beg2int128rc(const string& str){
	kmer res(0);
	for(int i(k-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[i]);
	}
	return res;
}


kmer graph3::end2int128rc(const string& str){
	kmer res(0);
	for(int i(k-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[str.size()-k+i]);
	}
	return res;
}


kmer graph3::end2int128(const string& str){
	kmer resEnd(0);
	for(uint i(0);i<k;++i){
		resEnd<<=2;
		resEnd+=chartoint(str[str.size()-k+i]);
	}
	return resEnd;
}


kmer graph3::rcb(kmer min){
	kmer resrcb(0);
	kmer offsetrcb(1);
	for(uint i(0); i<k;++i){
		resrcb+=(3-(min%4))<<(2*(k-1-i));
		min>>=2;
	}
	return resrcb;
}


void graph3::compaction(uint iL,  uint iR){
	if(iR!=iL){
		uint s1(unitigs[iL].size()),s2(unitigs[iR].size());
		bool b1(isNumber(unitigs[iL][0])),b2(isNumber(unitigs[iR][0]));
		if(b1 and b2){return compaction(stoi(unitigs[iL]),stoi(unitigs[iR]));}
		if(b1){return compaction(stoi(unitigs[iL]),iR);}
		if(b2){return compaction(iL,stoi(unitigs[iR]));}

		kmer beg1(beg2int128(unitigs[iL]));
		kmer end2(end2int128(unitigs[iR]));
		if(beg1==end2){
			unitigs[iR]+=(unitigs[iL].substr(k));
			unitigs[iL]=to_string(iR);
			return;
		}

		kmer endrc2(beg2int128rc(unitigs[iR]));
		if(beg1==endrc2){
			reverseinplace2(unitigs[iR]);
			unitigs[iR]+=(unitigs[iL].substr(k));
			unitigs[iL]=to_string(iR);
			return;
		}

		kmer beg2(rcb(endrc2));
		kmer end1(end2int128(unitigs[iL]));
		if(end1==beg2){
			unitigs[iL]+=(unitigs[iR].substr(k));
			unitigs[iR]=to_string(iL);
			return;
		}

		kmer begrc2(rcb(end2));
		if(end1==begrc2){
			unitigs[iL]+=(reverseinplace(unitigs[iR]).substr(k));
			unitigs[iR]=to_string(iL);
			return;
		}
	}
}


void graph3::debruijn(){
	sort(left.begin(),left.end(),comparator());
	sort(right.begin(),right.end(),comparator());
	uint iL(0),iR(0),sizeLeft(left.size()),sizeRight(right.size());
	left.push_back({0,(kmer)-1});
	right.push_back({0,(kmer)-1});
	kmerIndice kL,kR;
	while(iL!=sizeLeft and iR!=sizeRight){
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
				while(left[++iL].kmmer==kL.kmmer){}
			}else{
				while(right[++iR].kmmer==kR.kmmer){}
			}
		}
	}
}


bool graph3::output(uint i){return !isNumber(unitigs[i][0]);}


bool graph3::clear(){delete [] unitigs;return true;}


uint graph3::size(){return indiceUnitigs;};


void graph3::addtuple(tuple<string,uint,uint>& tuple){
	unitigs[indiceUnitigs]=move(get<0>(tuple));
	if(minimizer==(get<1>(tuple))){
		kmer kmer1(beg2int128(unitigs[indiceUnitigs]));
		kmer kmer2(rcb(kmer1));
		if(kmer1<kmer2){
			left.push_back(kmerIndice{indiceUnitigs,kmer1});
		}else{
			right.push_back(kmerIndice{indiceUnitigs,kmer2});
		}
	}
	if(minimizer==get<2>(tuple)){
		kmer kmer1(end2int128(unitigs[indiceUnitigs]));
		kmer kmer2(rcb(kmer1));
		if(kmer1<kmer2){
			right.push_back(kmerIndice{indiceUnitigs,kmer1});
		}else{
			left.push_back(kmerIndice{indiceUnitigs,kmer2});
		}
	}
	++indiceUnitigs;
}


// void compareUnitigs(const string& fileFa,const string& fileDot){
// 	uint a(0),b(0),c(0),d(0);
// 	unordered_set<string> setFa,setDot;
// 	ifstream streamFa(fileFa),streamDot(fileDot);
// 	string seq;
// 	getline(streamFa,seq);
// 	while (!streamFa.eof()) {
// 		getline(streamFa,seq,'>');
// 		seq=seq.substr(0,seq.size()-1);
// 		setFa.insert(seq);
// 		// cout<<seq<<endl;
// 		// cin.get();
// 		getline(streamFa,seq);
// 		++c;
// 	}
// 	cout<<1<<endl;
// 	while (!streamDot.eof()){
// 		getline(streamDot,seq);
// 		transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
// 		seq=seq.substr(0,seq.size()-1);
// 		setDot.insert(seq);
// 		// cout<<seq<<endl;
// 		// cin.get();
// 		++d;
// 	}
// 	cout<<2<<endl;
// 	for(auto it(setFa.begin());it!=setFa.end();++it){
// 		if(setDot.count(*it)==0){
// 			++a;
// 		}
// 	}
// 	cout<<3<<endl;
// 	for(auto it(setDot.begin());it!=setDot.end();++it){
// 		if(setFa.count(*it)==0){
// 			++a;
// 		}
// 	}
// 	cout<<a<<" "<<b<<endl;
// 	cout<<c<<" "<<d<<endl;
// }
//
//
// void compareKmers(const string& fileFa,const string& fileDot){
// 	uint k(31);
// 	string kmer;
// 	uint a(0),b(0),c(0),d(0);
// 	unordered_set<string> setFa,setDot;
// 	ifstream streamFa(fileFa),streamDot(fileDot);
// 	string seq,inter,nimp;
//
//
//
// 	// cout<<1<<endl;
// 	while (!streamFa.eof()) {
// 		getline(streamFa,nimp);
// 		// cout<<"nimp"<<nimp<<endl;
// 		getline(streamFa,seq);
// 		// cout<<"seq"<<seq<<endl;
// 		point:
// 		char c=streamFa.peek();
// 		if(c=='>'){
// 			point2:
// 			// seq=seq.substr(0,seq.size());
// 			// for(uint j(0);(j)<seq.size();++j){
// 			// 	if(seq[j]!='A' and seq[j]!='C' and seq[j]!='T' and seq[j]!='G'){
// 			// 		cout<<seq<<endl;
// 			// 		cout<<"lol"<<endl;
// 			// 		exit(0);
// 			// 	}
// 			// }
// 			for (uint i = 0; i+k <=seq.size(); ++i) {
// 				kmer=seq.substr(i,k);
// 				// cout<<kmer<<endl;
// 				kmer=getRepresent(kmer);
// 				// if(setDot.count(kmer)==0){
// 				// 	++a;
// 				// }
// 				setFa.insert(kmer);
// 			}
// 		}else{
// 			if(!streamFa.eof()){
// 				// cout<<"inter"<<endl;
// 				// cout<<seq<<endl;
// 				getline(streamFa,inter);
// 				// cout<<inter<<endl;
// 				seq+=inter;
// 				goto point;
// 			}else{
// 				// cout<<"lol2"<<endl;
// 				goto point2;
// 			}
// 		}
// 	}
// 	cout<<2<<endl;
//
// 	while (!streamDot.eof()){
// 		getline(streamDot,seq);
// 		seq=seq.substr(0,k);
// 		// cout<<seq<<endl;
// 		// cin.get();
// 		if(setFa.count(getRepresent(seq))==0){
// 			cout<<seq<<endl;
// 			++a;
// 		}
// 	}
//
// 	// while (!streamDot.eof()){
// 	// 	getline(streamDot,seq);
// 	// 	transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
// 	// 	seq=seq.substr(0,seq.size()-1);
// 	// 	// cout<<seq<<endl;
// 	// 	for (uint i = 0; i+k <=seq.size(); ++i) {
// 	// 		kmer=seq.substr(i,k);
// 	// 		// cout<<kmer<<endl;
// 	// 		kmer=getRepresent(kmer);
// 	// 		// setDot.insert(kmer);
// 	// 		if(setFa.count(kmer)==0){
// 	// 			++b;
// 	// 		}
// 	// 	}
// 	// 	// cout<<seq<<endl;
// 	// 	// cin.get();
// 	// 	// ++d;
// 	// }
// 	// for(auto it(setFa.begin());it!=setFa.end();++it){
// 	// 	if(setDot.count(*it)==0){
// 	// 		++a;
// 	// 	}
// 	// }
// 	cout<<3<<endl;
// 	// for(auto it(setDot.begin());it!=setDot.end();++it){
// 	// 	if(setFa.count(*it)==0){
// 	// 		++b;
// 	// 	}
// 	// }
// 	cout<<a<<" "<<b<<endl;
// 	cout<<c<<" "<<d<<endl;
// }
