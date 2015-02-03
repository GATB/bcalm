#include "ograph.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include "binSeq.h"

/*
 * constructs the compressed dBG from a list of arbitrary sequences
 */

using namespace std;

#ifdef SPARSE_HASH
    using namespace google;
#endif /* SPARSE_HASH */

static inline int nt2num(char c){
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

static inline int nt2num(char c, int pos){
	int offset = pos % 4;
	return (nt2num(c) + offset ) % 4;
}


static inline char num2nt(int num, int pos){
	int offset = pos % 4;
	return num2nt((num - offset + 4) % 4);
}

int shash(const string& s, int& previous_shash, unsigned int start_pos, int length){
    if (length == -1)
    {
        length = s.length();
    }
	int retval = 0, offset = 0;
    if (previous_shash != -1)
    {
        // shortcut, assuming consecutive m-mers
        retval = (previous_shash >> 2);
        offset = length - 1;
    }
	for (unsigned int pos = start_pos; pos+ offset < start_pos + length; pos++){
		//int num = nt2num(s[pos], pos - start_pos); // would need to be adapted to previous_hash
		int num = nt2num(s[pos + offset], 0);
		//retval = retval + pow(4, s.length() - pos - 1) * num;
        retval +=  (num  << (2*(pos + offset - start_pos))); // actually compute the reverse 2bit representation now
	}
    previous_shash = retval;

    // debug
    // printf("shash %s -> %0X\n",s.substr(start_pos,length).c_str(),retval);
	return retval;
}

string inverse_shash(int num, int len){
	string s(len, 'X');
    //~ int original_num = num;
	for (int pos = 0; pos < len; pos++){
		int val = num % 4;
		num = (num - val) / 4;
		s[pos] = num2nt(val);
	}

    // debug
    // printf("inverse shash %0X -> %s\n", original_num, s.c_str());
	return s;
}

uint32_t *m_mer_counts;

void init_m_mers_table(int m)
{
    m_mer_counts = new uint32_t[(int)pow(4,m)];
    //printf("Allocating %d MB for m-mers\n",(((int)pow(4,m)*sizeof(uint32_t))/1024)/1024);
    for (int i = 0; i < (int)pow(4,m); i++)
        m_mer_counts[i] = 0;
}

void count_m_mers(string str, int m, int k)
{
    int previous_shash = -1;
	for(int i(0) ;i < k-m+1 ; i+=1){
        m_mer_counts[shash(str, previous_shash, i, m)] ++;
    }
}

uint32_t *best_possible_hash_function;

void create_hash_function_from_m_mers(int m)
{
    int rg = pow(4,m);
    vector<pair<int, int> > counts;
    for (int i(0); i < rg; i++)
    {
        if (m_mer_counts[i] > 0)
            counts.push_back(make_pair(m_mer_counts[i],i));
    }

    printf("Sorting m-mer counts\n");
    sort(counts.begin(),counts.end());

    //~ cout<< counts[counts.size()-1].first << " " << inverse_shash(counts[counts.size()-1].second,m) <<endl;
    //~ cout<< counts[counts.size()-2].first << " " << inverse_shash(counts[counts.size()-2].second,m) <<endl;
    //~ cout<< counts[counts.size()-3].first << " " << inverse_shash(counts[counts.size()-3].second,m) <<endl;

    delete[] m_mer_counts; // fixme: restore this
    best_possible_hash_function = new uint32_t[rg];

    for (int i = 0; i < rg ; i++)
        best_possible_hash_function[i] = 0;

    for (unsigned int i = 0; i < counts.size(); i++)
        best_possible_hash_function[counts[i].second] =  i;
}

// not needed anymore
HashMap build_hash_map(int len){
	HashMap hm;
    return hm; // deactivate it

	for(int i(0);i<pow(4,len); i++){
		string s=inverse_shash(i,len);
		hm[s]=i;
	}
	hm.rehash(2 *(pow(4,len)));
	return hm;
}



int getash(const string& s, int& previous_shash, int start_pos = 0, int length = -1){
	//return ((*hm)[s]); // slow
    //return shash(s, start_pos, length);
    return best_possible_hash_function[shash(s, previous_shash, start_pos, length)];
}

int minimiserv(const string &node,const int &minimisersize){
    int previous_shash = -1;
	int minimiser_value(getash(node, previous_shash, 0, minimisersize)),vsub;
	for(uint64_t i(1);i<node.size()-minimisersize+1;i++){
		vsub=getash(node, previous_shash,i,minimisersize);
		if( minimiser_value > vsub){
			minimiser_value = vsub;
		}
	}
	//~ assert(minimiser_value<0);

    return(minimiser_value);
}

int minimiserrc(const string &node,const int &minimisersize){
	int h1, h2;
	h1 = minimiserv(node,minimisersize);
	h2 = minimiserv(reversecompletment(node),minimisersize);
    return (h1 > h2) ? h2 : h1;
}

int minimiserrc_openmp(const string &node,const int &minimisersize){
	string nodes[2];
	int resHash[2];

	nodes[0] = node;
	nodes[1] = reversecompletment(node);

	//~ #pragma omp parallel for
	for (int i=0; i<2; i++){
		resHash[i] = minimiserv(nodes[i],minimisersize);
	}

	if(resHash[0] < resHash[1]){
		return resHash[0];
	}
	else{
		return resHash[1];
	}
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

//read n character
string readn(ifstream *file,uint64_t n){
	string contents(n,'0');
	file->read(&contents[0],n);
	return(contents);
}



bool adjacent(const string& node1,const  string& node2,int k){
	return(node1.substr(node1.size()-k+1,k-1)==node2.substr(0,k-1));
}


static inline int chartoint(char c){
	switch(c){
		case 'A':
		return 0;
		case 'C':
		return 1;
		case 'G':
		return 2;
		case 'T':
		return 3;
		default:
		//~ cout<<"Problem with chartoint:"<<c<<endl;
		//~ assert(0);
		return 0;
	}
}

static inline char complement(char c){
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


void graph1::addvertex(const string& str){
	n++;
	unitigs.push_back(str);
	uint32_t i(unitigs.size()-1);
	uint64_t key(getkey(str));
	uint64_t keyrc(getkeyrevc(str));
	map.insert({key,i});
	maprev.insert({keyrc,i});
}

void graph1::addleftmin(int mini){
	leftmins.push_back(mini);
}

void graph1::addrightmin(int mini){
	rightmins.push_back(mini);
}

void graph1::debruijn(){
	neighbor=vector<neighbour> (n);
	string node,kmmer,kmmerr;
	uint64_t key,keyrc;
	for(uint64_t i(1);i<n;i++){
		node=unitigs[i];
		key=stringtoint(node.substr(node.size()-k+1,k-1));
		keyrc=stringtointc(node.substr(0,k-1));
		auto it(map.equal_range(key));
		for(auto j(it.first);j!=it.second;j++){
			//if k>32 collision can occur
			if(adjacent(node,unitigs[j->second],k)){
				neighbor[i].add(j->second,1);
				neighbor[j->second].add(i,4);
			}
		}
		it=(maprev.equal_range(key));
		for(auto j(it.first);j!=it.second;j++)
			if(adjacent(node,reversecompletment(unitigs[j->second]),k)){
				neighbor[i].add(j->second,2);
				neighbor[j->second].add(i,2);
			}
		it=(map.equal_range(keyrc));
		for(auto j(it.first);j!=it.second;j++)
			if(adjacent(reversecompletment(node),unitigs[j->second],k)){
				neighbor[i].add(j->second,3);
				neighbor[j->second].add(i,3);
			}
	}
	map.clear();
}


bool accordtomin(int min, int left_or_right_min){
	if(min == -1){
		return true;
	}

    if(left_or_right_min==min)
		return true;

	return false;

}

uint64_t graph1::becompacted(uint64_t nodeindice, int min, unsigned char *type){

	*type=0;
	string node=unitigs[nodeindice];

	if(node.empty())
		return 0;

    int leftmin = leftmins[nodeindice];
	int rightmin = rightmins[nodeindice];

	auto neigh(neighbor[nodeindice]);
	int one(neigh.nbtype(1)),two(neigh.nbtype(2)),three(neigh.nbtype(3)),four(neigh.nbtype(4));
	int in(three+four),out(one+two);
	if(out==1 && accordtomin(min,rightmin)){
		if(one==1){
			uint64_t sonindice(neigh.gtype(1));
			*type=1;
			if(neighbor[sonindice].nbtype(3)+neighbor[sonindice].nbtype(4)==1 && sonindice!=nodeindice)
				return sonindice;
		}
		else{
			uint64_t sonindice(neigh.gtype(2));
			*type=2;
			if(neighbor[sonindice].nbtype(2)+neighbor[sonindice].nbtype(1)==1 && sonindice!=nodeindice)
				return sonindice;
		}
	}
	if(in==1 && accordtomin(min,leftmin)){
		if(three==1){
			uint64_t sonindice(neigh.gtype(3));
			*type=3;
			if(neighbor[sonindice].nbtype(3)+neighbor[sonindice].nbtype(4)==1 && sonindice!=nodeindice)
				return sonindice;
		}
		else{
			uint64_t sonindice(neigh.gtype(4));
			*type=4;
			if(neighbor[sonindice].nbtype(1)+neighbor[sonindice].nbtype(2)==1 && sonindice!=nodeindice)
				return sonindice;
		}
	}
	return 0;
}

void graph1::reverse(int64_t with){
	string newnode(unitigs[with]);
	int64_t indice,type;
	if(newnode>(reversecompletment(newnode))){
		unitigs[with]=reversecompletment(newnode);

        int temp = leftmins[with];
        leftmins[with] = rightmins[with];
        rightmins[with] = temp;

		for(auto i(0);i<8;i++){
			indice=neighbor[with].list[i].first;
			type=neighbor[with].list[i].second;

			if(type==1){
				neighbor[indice].removep(with,4);
				neighbor[with].removep(indice,1);
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
				continue;
			}
			if(type==2){
				neighbor[indice].removep(with,2);
				neighbor[with].removep(indice,2);
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
				continue;
			}
			if(type==3){
				neighbor[indice].removep(with,3);
				neighbor[with].removep(indice,3);
				neighbor[with].add(indice,1);
				neighbor[indice].add(with,4);
				continue;
			}
			if(type==4){
				neighbor[indice].removep(with,1);
				neighbor[with].removep(indice,4);
				neighbor[with].add(indice,2);
				neighbor[indice].add(with,2);
				continue;
			}
		}
	}
}


void graph1::compact(uint64_t nodeindice,uint64_t with, unsigned char c){
	string newnode,node(unitigs[nodeindice]),son(unitigs[with]);
	uint64_t indice;
	if(nodeindice==with || node.empty() || son.empty())
		return;
	unsigned char type;
	switch(c){
		case 1:
		newnode=node+son.substr(k-1);
		unitigs[nodeindice]="";
		unitigs[with]=newnode;

        leftmins[with] = leftmins[nodeindice];
        //rightmins[with] = rightmins[with];

		for(int  i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			neighbor[indice].remove(nodeindice);
			if(indice==nodeindice){
				continue;
			}
			if(type==3 ){
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
			}
			if(type==4 ){
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
			}
		}
		break;

		case 2:
		newnode=node+reversecompletment(son).substr(k-1);
		unitigs[nodeindice]="";
		unitigs[with]=newnode;

        rightmins[with] = leftmins[with];
        leftmins[with] = leftmins[nodeindice];

		for(auto i(0);i<8;i++){
			indice=neighbor[with].list[i].first;
			type=neighbor[with].list[i].second;
			neighbor[indice].remove(with);
			neighbor[with].remove(indice);
			if(type==3){
				neighbor[with].add(indice,1);
				neighbor[indice].add(with,4);
			}
			if(type==4){
				neighbor[with].add(indice,2);
				neighbor[indice].add(with,2);
			}
		}
		for(auto i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			neighbor[indice].remove(nodeindice);
			if(type==3){
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
			}
			if(type==4 ){
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
			}
		}
		break;

		case 3:
		newnode=reversecompletment(node)+son.substr(k-1);
		unitigs[nodeindice]="";
		unitigs[with]=newnode;

        leftmins[with] = rightmins[nodeindice];
        //rightmins[with] = rightmins[with];


		neighbor[with].removep(nodeindice,3);
		for(auto i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			if(type==1){
				neighbor[indice].removep(nodeindice,4);
				neighbor[with].add(indice,3);
				neighbor[indice].add(with,3);
			}
			if(type==2 ){
				neighbor[indice].removep(nodeindice,2);
				neighbor[with].add(indice,4);
				neighbor[indice].add(with,1);
			}
		}
		break;

		case 4:
		newnode=son+node.substr(k-1);
		unitigs[nodeindice]="";
		unitigs[with]=newnode;

        // leftmins[with] = leftmins[with];
        rightmins[with] = rightmins[nodeindice];

		for(auto i(0);i<8;i++){
			indice=neighbor[nodeindice].list[i].first;
			type=neighbor[nodeindice].list[i].second;
			neighbor[indice].remove(nodeindice);
			if(type==1){
				neighbor[with].add(indice,1);
				neighbor[indice].add(with,4);
			}
			if(type==2 ){
				neighbor[with].add(indice,2);
				neighbor[indice].add(with,2);
			}
		}
		break;
	}
}


//Compact the graph but not the nodes that should be compacted in an other bucket
void graph1::compressh(int min){
	unsigned char type(0);
	uint64_t with;
	for(uint64_t nodeindice(1);nodeindice<n;nodeindice++){
		with=becompacted(nodeindice,min,&type);
		if(with!=0)
			compact(nodeindice,with,type);
	}
}

//Compact the graph but not the nodes that should be compacted in an other bucket
void graph1::compress(){
	HashMap hm;
	unsigned char type(0);
	uint64_t with;
	for(uint64_t nodeindice(1);nodeindice<n;nodeindice++){
		with=becompacted(nodeindice,-1,&type);
		if(with!=0)
			compact(nodeindice,with,type);
	}
}



//import graph from file
void graph1::importg(const char *name){
	cout<<"importg"<<endl;
	ifstream fichier(name,ios::in);
	string line,vertex;
	if(fichier)
		while(!fichier.eof()){
			getline(fichier,line);
			int64_t lp(0);
			for(unsigned int i(0);i<line.size();i++){
				if(line[i]==';'){
					addvertex(line.substr(lp,i-lp));
					lp=i+1;
				}
			}
		}
	else
		cerr<<"no file"<<endl;
}



int graph1::weight(){
	int res(0);
	for(auto k(unitigs.begin());k!=unitigs.end();k++)
		res+=k->size();
	return res;
}



void graph1::print(const char *name){
	ofstream fichier(name, ios::out | ios::trunc);
	if(fichier)
		for(uint64_t i(1);i<n;i++){
			string s(unitigs[i]);
			if(s!="")
				fichier<<s<<";"<<endl;
		}
}

void graph1::printedges(const char *name){
	ofstream fichier(name, ios::out | ios::trunc);
	if(fichier){
		fichier << "digraph test {" <<endl;
		for(uint64_t i(1);i<n;i++){
			string s(unitigs[i]);
			if(s!=""){
				fichier<<s<<";"<<endl;
				//~ auto v=neighbor[i].son;
				//~ for(auto j=v.begin();j!=v.end();j++)
					//~ if(j->first!=0)
						//~ if(nodes[j->first]!="")
							//~ fichier<<s<<"->"<<nodes[j->first]<<";"<<endl;
			}
		}
		fichier << "}"<<endl;
		fichier.close();
	}
}


uint64_t graph1::getkey(const string& str){
	return stringtoint(str.substr(0,k-1));
}


uint64_t graph1::getkeyrevc(const string& str){
	return stringtointc(str.substr(str.size()-k+1,k-1));
}


void neighbour::add(uint64_t p,unsigned char c){
	for(int i(0);i<8;i++){
		if(list[i].first==0  ){
			list[i]=make_pair(p,c);
			return;
		}
		if(list[i].first==p && list[i].second==c )
			return;
	}
}



uint64_t neighbour::nbtype(unsigned char c){
	uint64_t ret(0);
	for(int i(0);i<8;i++)
		if(list[i].second==c)
			ret++;
	return ret;
}




uint64_t neighbour::gtype(unsigned char c){
	for(int i(0);i<8;i++)
		if(list[i].second==c)
			return list[i].first;
	cout<<"Bug with neighbour"<<endl;
	return 0;
}




unsigned char neighbour::remove(uint64_t v){
	for(int i(0);i<8;i++)
		if(list[i].first==v){
			list[i].first=0;
			list[i].second=0;
			return 0;
		}
	return 0;
}

unsigned char neighbour::removep(uint64_t v,unsigned char c){
	for(int i(0);i<8;i++)
		if(list[i].first==v && list[i].second==c){
			list[i].first=0;
			list[i].second=0;
			return 0;
		}
	return 0;
}


unsigned char neighbour::removetype(unsigned char c){
	for(int i(0);i<8;i++)
		if(list[i].second==c)
			list[i].first=0;
	return 0;
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


void graph2::addvertex(const string& unitig){
	unitigs.push_back(unitig);
	size_t i=unitigs.size()-1;
	if(leftmins[i]){
		string beg(unitigs[i].substr(0,k));
		uint64_t leftKmer1(stringtoint(beg));
		uint64_t leftKmer2(stringtointc(beg));
		if(leftKmer1<leftKmer2){
			if(left2unitig.count(leftKmer1)>0){
				left2unitig[leftKmer1]=0;
			}else{
				left2unitig[leftKmer1]=i;
			}
		}else{
			if(right2unitig.count(leftKmer2)>0){
				right2unitig[leftKmer2]=0;
			}else{
				right2unitig[leftKmer2]=i;
			}
		}
	}

	if(rightmins[i]){
		string end(unitigs[i].substr(unitigs[i].size()-k,k));
		uint64_t rightKmer1(stringtoint(end));
		uint64_t rightKmer2(stringtointc(end));
		if(rightKmer1<rightKmer2){
			if(right2unitig.count(rightKmer1)>0){
				right2unitig[rightKmer1]=0;
			}else{
				right2unitig[rightKmer1]=i;
			}
		}else{
			if(left2unitig.count(rightKmer2)>0){
				left2unitig[rightKmer2]=0;
			}else{
				left2unitig[rightKmer2]=i;
			}
		}
	}
}

bool kmerIndiceCompare(kmerIndice a ,kmerIndice b) { return a.kmmer < b.kmmer; }

bool isNumber(const string& str){
	switch (str[0]){
	case 'a': return false;
		break;
	case 'c': return false;
		break;
	case 'g': return false;
		break;
	case 't': return false;
		break;
	case 'A': return false;
		break;
	case 'C': return false;
		break;
	case 'G': return false;
		break;
	case 'T': return false;
		break;
	default: return true;
}
}



void graph3::compaction(uint32_t iL, uint32_t iR){
	string seq1(unitigs[iL]),seq2(unitigs[iR]);
    size_t s1(seq1.size()),s2(seq2.size());

    bool b1(isNumber(seq1)),b2(isNumber(seq2));
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
		unitigs[iL]=seq2+seq1.substr(k);
		unitigs[iR]=to_string(iL);
        return;
    }else{
        if(rc2.substr(s2-k,k)==seq1.substr(0,k)){
		unitigs[iL]=rc2+seq1.substr(k);
		unitigs[iR]=to_string(iL);
		return;
        }
    }

    if(seq2.substr(0,k)==seq1.substr(s1-k,k)){
        unitigs[iL]=seq1+seq2.substr(k);
		unitigs[iR]=to_string(iL);
		return;
    }else{
        if(rc1.substr(s1-k,k)==seq2.substr(0,k)){
			unitigs[iL]=rc1+seq2.substr(k);
			unitigs[iR]=to_string(iL);
			return;
        }
    }

    if(rc1.substr(0,k)==seq2.substr(s2-k,k)){
            unitigs[iL]=seq2+rc1.substr(k);
			unitigs[iR]=to_string(iL);
			return;
    }else{
        if(rc2.substr(s2-k,k)==rc1.substr(0,k)){
            unitigs[iL]=rc2+rc1.substr(k);
			unitigs[iR]=to_string(iL);
			return;
        }
    }

    if(rc2.substr(0,k)==seq1.substr(s1-k,k)){
        unitigs[iL]=seq1+rc2.substr(k);
		unitigs[iR]=to_string(iL);
		return;
    }else{
        if(rc1.substr(s1-k,k)==rc2.substr(0,k)){
            unitigs[iL]=rc1+rc2.substr(k);
			unitigs[iR]=to_string(iL);
			return;
        }
    }
    cout<<"WUT"<<endl;
    cout<<"+"<<seq1<<"+ +"<<seq2<<"+"<<endl;
    cin.get();
}



void graph3::debruijn(){
        sort(left.begin(),left.end(),kmerIndiceCompare);
        sort(right.begin(),right.end(),kmerIndiceCompare);

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
                                   kmer2Indice k2i;
                                k2i.indiceL=kL.indice;
                                k2i.indiceR=kR.indice;
                                compactions.push_back(k2i);
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
    left.clear();
        right.clear();
}



void graph3::compress(){
	for(size_t i(0);i<compactions.size();++i){
		kmer2Indice k2i(compactions[i]);
		uint32_t iL(k2i.indiceL);
		uint32_t iR(k2i.indiceR);
		compaction(iL,iR);
	}

	for(size_t i(0);i<unitigs.size();++i){
		if(isNumber(unitigs[i])){
			unitigs[i]="";
		}
	}
}

void graph3::addvertex(const string& unitig){
	unitigs.push_back(unitig);
	size_t i=unitigs.size()-1;
	if(leftmins[i]){
		string beg(unitigs[i].substr(0,k));
		__uint128_t leftKmer1(stringtoint128(beg));
		__uint128_t leftKmer2(stringtointc128(beg));
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
		string end(unitigs[i].substr(unitigs[i].size()-k,k));
		__uint128_t rightKmer1(stringtoint128(end));
		__uint128_t rightKmer2(stringtointc128(end));
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

void graph3::addleftmin(int min){
	if(min==minimizer){
		leftmins.push_back(true);
	}else{
		leftmins.push_back(false);
	}
}

void graph3::addrightmin(int min){
	if(min==minimizer){
		rightmins.push_back(true);
	}else{
		rightmins.push_back(false);
	}
}

void graph2::addleftmin(int min){
	if(min==minimizer){
		leftmins.push_back(true);
	}else{
		leftmins.push_back(false);
	}
}

void graph2::addrightmin(int min){
	if(min==minimizer){
		rightmins.push_back(true);
	}else{
		rightmins.push_back(false);
	}
}

uint32_t graph2::leftUnique(uint64_t i){
	if(left2unitig.count(i)>0){
		//~ return(left2unitig.at(i));
		return(left2unitig[i]);
	}
	return 0;
}

uint32_t graph2::rightUnique(uint64_t i){
	if(right2unitig.count(i)>0){
		//~ return(right2unitig.at(i));
		return(right2unitig[i]);
	}
	return 0;
}

uint32_t graph2::goEnd(uint32_t i){
	//~ if(usedUnitigs.count(i)==0){
		if(rightmins[i]){
			if(unitigs[i].size()!=0){
				string end(unitigs[i].substr(unitigs[i].size()-k,k));
				uint64_t rightKmer1(stringtoint(end));
				uint64_t rightKmer2(stringtointc(end));
				if(rightKmer1<rightKmer2){
					if(rightUnique(rightKmer1)!=0){
						//~ cout<<"rightunique"<<endl;
						return(leftUnique(rightKmer1));
					}
				}else{
					if(leftUnique(rightKmer2)!=0){
						//~ cout<<"leftunique"<<endl;
						return(rightUnique(rightKmer2));
					}
				}
			}
		}
	return 0;
}

uint32_t graph2::goBeg(uint32_t i){
	//~ if(usedUnitigs.count(i)==0){
		if(leftmins[i]){
			if(unitigs[i].size()!=0){
				string beg(unitigs[i].substr(0,k));
				uint64_t leftKmer1(stringtoint(beg));
				uint64_t leftKmer2(stringtointc(beg));
				if(leftKmer1<leftKmer2){
					if(leftUnique(leftKmer1)!=0){
						//~ cout<<"leftunique"<<endl;
						return(rightUnique(leftKmer1));
					}
				}else{
					if(rightUnique(leftKmer2)!=0){
						//~ cout<<"lrightunique"<<endl;
						return(leftUnique(leftKmer2));
					}
				}
			}
		}
	//~ }
	return 0;
}

void graph2::debruijn(){
    return; // we don't do anything in here anymore
	for(uint32_t i=1;i<unitigs.size();++i){

	}
}

uint32_t graph1::size(){
	return unitigs.size();
}

uint32_t graph2::size(){
	return unitigs.size();
}

uint32_t graph3::size(){
	return unitigs.size();
}

uint32_t graph4::size(){
	return unitigs.size();
}

void graph2::compress(){
	for(uint32_t i=1;i<unitigs.size();++i){
		if(unitigs[i].size()!=0){
			uint32_t next(goEnd(i));
			if(next!=0){
				if(goBeg(i)==0){
					chainCompaction2(i,unitigs[i],next);
				}
			}
			else{
				next=goBeg(i);
				if(next!=0){
					chainCompaction2(i,unitigs[i],next);
				}
			}
		}
	}
}

void graph2::chainCompaction2(uint32_t i, string unitig, uint32_t next){
	unordered_set<uint32_t> LocalusedUnitigs;
	LocalusedUnitigs.insert(i);
	uint32_t last(next);
	string chain(compaction(unitig,unitigs[next],k));
	bool change(false);
	while(chain!=unitig){
		change=true;
		unitig=chain;
		last=next;
		LocalusedUnitigs.insert(next);
		next=goEnd(next);
		if(next!=0 and LocalusedUnitigs.count(next)==0){
			chain=compaction2(unitig,unitigs[next],k);
		}else{
			next=last;
			next=goBeg(next);
			if(next!=0 and LocalusedUnitigs.count(next)==0){
				chain=compaction2(unitig,unitigs[next],k);
			}
		}
	}
	if(change){
		//~ cout<<"no"<<endl;
		//~ if(getRepresent(chain.substr(0,k))==getRepresent(unitigs[last].substr(0,k))){
			//~ leftmins[i]=leftmins[last];
		//~ }
		//~ if(getRepresent(chain.substr(chain.size()-k,k))==getRepresent(unitigs[last].substr(unitigs[last].size()-k,k))){
			//~ rightmins[i]=rightmins[last];
		//~ }
		//~ if(getRepresent(chain.substr(0,k))==getRepresent(unitigs[last].substr(unitigs[last].size()-k,k))){
			//~ leftmins[i]=rightmins[last];
		//~ }
		//~ if(getRepresent(chain.substr(chain.size()-k,k))==getRepresent(unitigs[last].substr(0,k))){
			//~ rightmins[i]=leftmins[last];
		//~ }
		unitigs[i]=chain;
		LocalusedUnitigs.erase(i);
		for (const std::uint32_t& x: LocalusedUnitigs){
			unitigs[x]="";
		}
	}
}







void graph4::debruijn(){
        sort(left.begin(),left.end(),kmerIndiceCompare);
        sort(right.begin(),right.end(),kmerIndiceCompare);

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
                                   kmer2Indice k2i;
                                k2i.indiceL=kL.indice;
                                k2i.indiceR=kR.indice;
                                compactions.push_back(k2i);
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
    left.clear();
    right.clear();
}



void graph4::compress(){
	for(size_t i(0);i<compactions.size();++i){
		kmer2Indice k2i(compactions[i]);
		uint32_t iL(k2i.indiceL);
		uint32_t iR(k2i.indiceR);
		//~ cout<<iL<<" "<<iR<<endl;
		compaction(iL,iR);
	}

	for(size_t i(0);i<unitigs.size();++i){
		if(unitigs[i].isNumber){
			unitigs[i].clear();
		}
	}
}

void graph4::addvertex(const string& unitigstr){
	//~ cout<<"ADD "<<endl<<unitigstr<<endl;
	//~ cout<<"------------------------------------------------------------------------------------------------------------"<<endl;
	binSeq unitig(unitigstr);
	//~ cout<<unitig.size()<<endl;
	unitigs.push_back(unitig);
	uint32_t i(unitigs.size()-1);

	if(leftmins[i]){
		binSeq beg(unitig.getBegin(k));
		//~ cout<<beg.str()<<endl;
		uint64_t leftKmer1(beg.getInt());
		//~ cout<<leftKmer1<<endl;
		beg.reverse();
		//~ cout<<beg.str()<<endl;
		uint64_t leftKmer2(beg.getInt());
		//~ cout<<leftKmer2<<endl;
		kmerIndice ki;
		ki.indice=i;
		if(leftKmer1<leftKmer2){
			ki.kmmer=leftKmer1;
			//~ cout<<"left"<<leftKmer1<<endl;
			left.push_back(ki);
		}else{
			ki.kmmer=leftKmer2;
			//~ cout<<"right"<<leftKmer2<<endl;
			right.push_back(ki);
		}
	}

	if(rightmins[i]){
		binSeq end(unitig.getEnd(k));
		//~ cout<<end.str()<<endl;
		uint64_t rightKmer1(end.getInt());
		end.reverse();
		//~ cout<<rightKmer1<<endl;
		//~ cout<<end.str()<<endl;
		uint64_t rightKmer2(end.getInt());
		//~ cout<<rightKmer2<<endl;
		kmerIndice ki;
		ki.indice=i;
		if(rightKmer1<rightKmer2){
			ki.kmmer=rightKmer1;
			//~ cout<<"right"<<rightKmer1<<endl;
			right.push_back(ki);
		}else{
			ki.kmmer=rightKmer2;
			//~ cout<<"left"<<rightKmer2<<endl;
			left.push_back(ki);
		}
	}

	//~ cout<<"ADDED"<<endl<<endl;;
}

void graph4::addleftmin(int min){
	if(min==minimizer){
		leftmins.push_back(true);
	}else{
		leftmins.push_back(false);
	}
}

void graph4::addrightmin(int min){
	if(min==minimizer){
		rightmins.push_back(true);
	}else{
		rightmins.push_back(false);
	}
}



void graph4::compaction(uint32_t iL, uint32_t iR){

    bool b1(unitigs[iL].isNumber),b2(unitigs[iR].isNumber);
    if(b1 && b2){
		compaction(unitigs[iL].getNumber(),unitigs[iR].getNumber());
		return;
	}
    if(b1){
		compaction(unitigs[iL].getNumber(),iR);
		return;
	}
	if(b2){
		compaction(iL,unitigs[iR].getNumber());
		return;
	}

	if(unitigs[iL].size()<unitigs[iR].size()){
		swap(iR,iL);
	}

	//~ binSeq seq1(unitigs[iL]),seq2(unitigs[iR]);

		binSeq end1(unitigs[iL].getEnd(k));
		binSeq beg2(unitigs[iR].getBegin(k));
		if(end1.getInt()==beg2.getInt()){
			unitigs[iL].add(unitigs[iR].sub(k));
			unitigs[iR]=binSeq(iL);
			return;
		}

		binSeq end2(unitigs[iR].getEnd(k));
		if(end1.getInt()==end2.getReverse().getInt()){
			unitigs[iL].add(unitigs[iR].getReverse().sub(k));
			unitigs[iR]=binSeq(iL);
			return;
		}

		binSeq beg1(unitigs[iL].getBegin(k));
		if(beg1.getInt()==end2.getInt()){
			unitigs[iR].add(unitigs[iL].sub(k));
			unitigs[iL]=binSeq(iR);
			return;
		}

		if(beg1.getInt()==beg2.getReverse().getInt()){
			unitigs[iR].reverse();
			unitigs[iR].add(unitigs[iL].sub(k));
			unitigs[iL]=binSeq(iR);
			return;
		}


    cout<<"WUT"<<endl;
    cout<<"+"<<unitigs[iL].str()<<"+ +"<<unitigs[iR].str()<<"+"<<endl;
}






















