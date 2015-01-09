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

using namespace std;

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



string reversecomplement(const string& dna) // code taken from TakeABreak
{
	string revComp= "";

	for (string::const_reverse_iterator it = dna.rbegin(); it != dna.rend(); it++) {
		switch(*it) {
			case 'a' :                revComp += "t";                break;
			case 't' :                revComp += "a";                break;
			case 'c' :                revComp += "g";                break;
			case 'g' :                revComp += "c";                break;
			case 'A' :                revComp += "T";                break;
			case 'T' :                revComp += "A";                break;
			case 'C' :                revComp += "G";                break;
			case 'G' :                revComp += "C";                break;
		}
	}
    return revComp;
}

/* nifty function to highlight a substring for printing to terminal */
string debug_highlight(string s, string motif)
{
    size_t pos = s.find(motif);
    if (pos == string::npos)
        pos = s.find(reversecomplement(motif));
    if (pos == string::npos)
    {
        return s;
    }

    return s.substr(0,pos) + "\033[1;31m" + s.substr(pos, motif.size()) + "\033[0m" + s.substr(pos + motif.size(), s.size() - pos - motif.size()); 
}


void rcGlueEntry(GlueEntry & entry) {
	int k = entry.lkmer.size();
	entry.seq = reversecomplement(entry.seq);
	std::swap(entry.lmark, entry.rmark);
	entry.lkmer = entry.seq.substr(0,k);
	entry.rkmer = entry.seq.substr(entry.seq.length() - k, k);
}

GlueEntry GlueEntryCompactNaive::getEntry() {
	GlueEntry e;
	//if (raw.length() < kmerSize + 2) { cout << "problem: " << raw << ".\n"; exit(1); }
	if (raw.length() < kmerSize + 2) { //an empty entry
		return e;
	}
	e.seq = raw.substr(0, raw.length() - 2);
	e.lkmer = e.seq.substr(0, kmerSize);
	e.rkmer = e.seq.substr(e.seq.length() - kmerSize, kmerSize);
	e.lmark = raw[raw.length() - 2];
	e.rmark = raw[raw.length() - 1];
	return e;
}


bool GlueStorage::find (string key, GlueEntry e1, GlueEntry e2) {
	findIt  = glueMap.find (key);
	if ( findIt == glueMap.end() ) return false;
	GlueEntryCompactNaive e1c (findIt->second.first, kmerSize);
	GlueEntryCompactNaive e2c (findIt->second.second, kmerSize);
	e1 = e1c.getEntry();
	e2 = e2c.getEntry();
	return true;
}

bool GlueStorage::getFirst(GlueEntry & e1, GlueEntry & e2) {
	curIt = glueMap.begin();
	return derefIt(curIt, e1, e2);

}
bool GlueStorage::getNext(GlueEntry & e1, GlueEntry & e2) {
	curIt++;
	return derefIt(curIt, e1,e2);
}

void GlueStorage::insertAtKey(string key, GlueEntry e1, GlueEntry e2) {
	string s1 = GlueEntryCompactNaive(e1, kmerSize).getRaw();
	string s2 = GlueEntryCompactNaive(e2, kmerSize).getRaw();
	glueMap[key] = make_pair(s1, s2);
}


void GlueStorage::insertAfterFind(GlueEntry e1, GlueEntry e2) {
	if (findIt == glueMap.end()) {
		cout << "GlueStorage::insertAfterFind has findIt pointing to the end of glueMap...SHOULD not BE.\n";
		exit(1);
	}
	string s1 = GlueEntryCompactNaive(e1, kmerSize).getRaw();
	string s2 = GlueEntryCompactNaive(e2, kmerSize).getRaw();
	findIt->second.first = s1;
	findIt->second.second = s2;
}

void GlueStorage::insertAtCurIt(GlueEntry e1, GlueEntry e2) {
	if (curIt == glueMap.end()) {
		cout << "GlueStorage::insertAtCurIt has curIt pointing to the end of glueMap...SHOULD not BE.\n";
		exit(1);
	}
	string s1 = GlueEntryCompactNaive(e1, kmerSize).getRaw();
	string s2 = GlueEntryCompactNaive(e2, kmerSize).getRaw();
	curIt->second.first = s1;
	curIt->second.second = s2;
}


void GlueStorage::cleanup() {
	GlueEntry e1, e2;
	for (auto it = glueMap.begin(); it != glueMap.end(); ) {
		derefIt(it, e1, e2);
		if ((e1.seq == "") && (e2.seq == "")) {
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

void GlueStorage::updateMemStats() {
	maxEntries = std::max(maxEntries, glueMap.size());
	totEntries += glueMap.size();

	size_t size = 0;
	for (auto it_gS = glueMap.begin(); it_gS != glueMap.end(); it_gS++)
	{
		size += sizeof(it_gS->first) + sizeof(it_gS->second.first) + sizeof(it_gS->second.second);
		size += it_gS->first.size() + it_gS->second.first.size() + it_gS->second.second.size();
	}

	maxSize = std::max(maxSize, size);
	totSize += size;
	numDataPoints++;
}

bool GlueStorage::derefIt (GlueMap::const_iterator it, GlueEntry e1, GlueEntry e2) {
	if (it == glueMap.end()) return false;
	e1 = GlueEntryCompactNaive(it->second.first, kmerSize).getEntry();
	e2 = GlueEntryCompactNaive(it->second.second, kmerSize).getEntry();
	return true;
}


void GlueStorage::printMemStats() {
	if (numDataPoints == 0) {
		cout << "GlueStorage: no data points to output memory stats.\n";
	} else {
		cout << "GlueStorage memory stats: max entries: " << add_commas(maxEntries) << " avg entries: " << add_commas(totEntries/numDataPoints) << " max size: " << add_commas(maxSize) 
			<< "b avg size: " << add_commas(totSize / numDataPoints) << "b\n";
	}
}


void Glue::insert(GlueEntry newEntry) { 
	string key = newEntry.rkmer;

	if (key.compare(reversecomplement(key)) <= 0) {
		rcGlueEntry(newEntry);
		key = reversecomplement(key);
	}

	GlueEntry e1, e2;
	if (!glueStorage.find(key, e1, e2)) {
		glueStorage.insertAtKey(key, newEntry, GlueEntry());
	} else {
		if ((!e1.lmark && !e1.rmark) || (e1.seq == "")) {
			glueStorage.insertAfterFind(newEntry, GlueEntry());
		} else {
			if (e2.seq == "") {
				glueStorage.insertAfterFind(e1, newEntry);
			} else {
				printf("huh? kmer %s (%s,%s) inserting %s\n",key.c_str(), debug_highlight(e1.seq, key).c_str(), debug_highlight(e2.seq, key).c_str(), debug_highlight(newEntry.seq, key).c_str());
				exit(1);
			}
		}
	}

}


void Glue::output(string seq)
{
	Sequence s (Data::ASCII);
	s.getData().setRef ((char*)seq.c_str(), seq.size());
	out.insert(s);
}

void Glue::glue()
{
	size_t k = kmerSize;

	
	GlueEntry query, match;
	for (bool active = glueStorage.getFirst(query, match); active; active = glueStorage.getNext(query, match)) {
		if ((query.seq == "") || (match.seq == "")) 
			continue; // not yet ready to process this
		if ((!query.lmark && !query.rmark)) {
			if (!match.lmark && !match.rmark) {
				printf ("Assert failed!\n"); 
				exit(1);
			}
			continue; //already output
		}
		if (!match.lmark  || (query.rkmer.compare(match.lkmer) != 0)) {
			std::swap(query, match);
		}
		if (!match.lmark) {
			continue; // nothing to do here
		}
		if (query.rkmer.compare(match.lkmer) == 0) {
			glueStorage.insertAtCurIt(GlueEntry(), GlueEntry());  //effectively removes the entry from the table
			string gluedStr = query.seq + match.seq.substr(k, match.seq.size() - k);
			if (!query.lmark && !match.rmark) {
				output(gluedStr);
			} else {
				GlueEntry e;
				e.seq = gluedStr;
				e.lmark = query.lmark;
				e.rmark = query.rmark;
				e.lkmer = e.seq.substr(0,k);
				e.rkmer = e.seq.substr(e.seq.length() - k, k);
				insert(e);
			}
		} else {
			cout << "uh oh, not matching" << endl;
			cout << "  query: " << debug_highlight(query.seq, query.rkmer) <<  "\n";
			cout << "  match: " << debug_highlight(match.seq, match.lkmer) << std::endl;
		}
	}

	glueStorage.cleanup();
}
/* 
 * class to support compact representation of a glueentry, as should be stored in memory
 * this version will support a binary representation
 * NOT ACTIVE YET
 *
class GlueEntryCompactBinary{
	public:
		const static int MAX_SIZE = 1000;
		bool raw[MAX_SIZE];
		int rawLen;

		GlueEntryCompact(bool * _raw, int _rawLen) : rawLen(_rawLen) { 
			for (int i = 0; i < rawLen; i++) {
				raw[i] = _raw[i];
			}
		}

		GlueEntryCompact(string seq, bool leftmark, bool rightmark) {
			rawLen = seq.length() * 2 + 2;
			for (int i = 0; i < seq.length(); i++) {
				if ((seq[i] == 'A') || (seq[i] == 'C')) {
					raw[2*i] = 0;
				} else {
					raw[2*i] = 1;
				}
				if ((seq[i] == 'A') || (seq[i] == 'G')) {
					raw[2*i + 1] = 0;
				} else {
					raw[2*i + 1] = 1;
				}
			}
			raw[seq.length() * 2] = leftmark;
			raw[seq.length() * 2 + 1] = rightmark;
		}


		char bits2char (bool leftbit, bool rightbit) {
			if (leftbit == 0 && rightbit == 0) return 'A';
			else if (leftbit == 0 && rightbit == 1) return 'C';
			else if (leftbit == 1 && rightbit == 0) return 'G';
			else return 'T';
		}


		GlueEntry getEntry() {
			GlueEntry entry;
			for (int  i = 0; i < rawLen - 2; i += 2) {
				entry.seq.push_back(bits2char(raw[i], raw[i+1]));
			}
			entry.lkmer = entry.seq.substr(0,k);
			entry.rkmer = entry.seq.substr(entry.seq.length() - k, k);
			entry.lmark = raw[rawLen - 2];
			entry.rmark = raw[rawLen - 1];
		}
};
*/



/*


   OLD GLUE CODE
   OLD GLUE CODE
   OLD GLUE CODE
   OLD GLUE CODE

   string int2string(int c)
   {
   if (c == 0)
   return "0";
   else if (c == 1)
   return "1";
   else
   return "2";
   }



   void Glue::remove(string seq)
   {
// extremely naive for now
string leftkmer = seq.substr(0,kmerSize);
string rightkmer = seq.substr(seq.size() - kmerSize - 2, kmerSize);

string min = std::min(leftkmer, reversecomplement(leftkmer));
if (seq_only(glueStrings[min].first) == seq_only(seq) \
|| seq_only(glueStrings[min].first) == reversecomplement(seq_only(seq)))
{
glueStrings[min].first = glueStrings[min].second;
glueStrings[min].second = "";
}
if (seq_only(glueStrings[min].second) == seq_only(seq) \
|| seq_only(glueStrings[min].second) == reversecomplement(seq_only(seq)))
{
glueStrings[min].second = "";
}

min = std::min(rightkmer, reversecomplement(rightkmer));
if (seq_only(glueStrings[min].first) == seq_only(seq) \
|| seq_only(glueStrings[min].first) == reversecomplement(seq_only(seq)))
{
glueStrings[min].first = glueStrings[min].second;
glueStrings[min].second = "";
}
if (seq_only(glueStrings[min].second) == seq_only(seq) \
|| seq_only(glueStrings[min].second) == reversecomplement(seq_only(seq)))
{
glueStrings[min].second = "";
}

}




void Glue::now_really_addhash(string key, string value)
{
GlueMap::const_iterator got = glueStrings.find (key);
if ( got == glueStrings.end() )
{
glueStrings[key].first = value;
glueStrings[key].second = "";
}
else
{
if (glueStrings[key].first == "-1" || glueStrings[key].first == "")
{
glueStrings[key].first = value;
glueStrings[key].second = "";
}
else
{
	if  (glueStrings[key].second != "")
	{
		printf("huh? kmer %s (%s,%s) inserting %s\n",key.c_str(), debug_highlight(glueStrings[key].first,key).c_str(), debug_highlight(glueStrings[key].second,key).c_str(), debug_highlight(value,key).c_str());
		exit(1);
	}
	glueStrings[key].second = value;
}
}

}
*/

/*
   void Glue::aux_addhash(string seq, string kmer, bool leftmark, bool rightmark, bool left)
   {
   if (kmer.compare(reversecomplement(kmer)) <= 0)
   now_really_addhash(kmer, seq + int2string(leftmark) + int2string(rightmark));
   else
   now_really_addhash(reversecomplement(kmer), reversecomplement(seq) + int2string(rightmark) + int2string(leftmark));
   return;

// this was the former strategy
if (kmer.compare(reversecomplement(kmer)) <= 0)
{
if (left)
glueStrings[kmer].second = seq + int2string(leftmark) + int2string(rightmark);
else
glueStrings[kmer].first  = seq + int2string(leftmark) + int2string(rightmark);
}
else
{
if (left)
glueStrings[reversecomplement(kmer)].first = reversecomplement(seq) + int2string(rightmark) + int2string(leftmark);
else
glueStrings[reversecomplement(kmer)].second = reversecomplement(seq) + int2string(rightmark) + int2string(leftmark);
}
}

void Glue::insert(string seq, bool leftmark, bool rightmark)
{
//if (leftmark)
{
// don't yet know how to index the hash table with kmers
//ModelCanon::Kmer kmer = model.codeSeed(seq.substr(0,kmerSize).c_str(), Data::ASCII); 
string kmer = seq.substr(0,kmerSize);
aux_addhash(seq, kmer, leftmark, rightmark, true);
}
//if (rightmark)
if (seq.size() > kmerSize)
{
//ModelCanon::Kmer kmer = model.codeSeed(seq.substr(seq.size() - kmerSize,kmerSize).c_str(), Data::ASCII); 
string kmer = seq.substr(seq.size() - kmerSize, kmerSize);
aux_addhash(seq, kmer, leftmark, rightmark, false);
}
}

void Glue::glue()
{
size_t k = kmerSize;

// mostly copypaste of code by Colleen
for (auto it_gS = glueStrings.begin(); it_gS != glueStrings.end(); it_gS++)
{
string query = (it_gS->second).first;   // "query" is always the first node, should be left side of final node
string match = (it_gS->second).second;

if (query == "" || match == "")
continue; // not yet ready to process this

assert((query.compare("-1") == 0 && match.compare("-1") == 0) || (query.compare("-1") != 0 && match.compare("-1") != 0)); // should never have just one string of pair -1
// fixme: for some reason assert does not work in this project 
if (query != "-1" && match=="-1") {printf("assert failed!\n"); exit(1);}

// already output
if (query == "-1")
{
continue;
}

string queryNode = query.substr(0, query.size()-2); // remove leftmark and rightmark
string rightKmer = queryNode.substr(queryNode.size() - k, k);                       //last k characters of the query node
string leftKmer = queryNode.substr(0, k);                                           //first k characters of the query node
string queryLeftMark = query.substr(query.size() - 2, 1);
string queryRightMark = query.substr(query.size() - 1, 1);

if (debug)
	cout << "query is: " << queryNode << queryLeftMark << queryRightMark << "\n";
if (debug)
	cout << "rightKmer is: " << rightKmer << "\n";

	string matchNode = match.substr(0, match.size()-2);
if (debug)
	cout << "matchnode: '" <<  matchNode << "' match: " << match << std::endl;

	string matchLeftKmer = matchNode.substr(0, k);
	string matchRightKmer = matchNode.substr(matchNode.size() - k, k);
	string matchLeftMark = match.substr(match.size() - 2, 1);
	string matchRightMark = match.substr(match.size() - 1, 1);

	if (matchLeftMark != "1" || rightKmer.compare(matchLeftKmer) != 0) // try swapping
{
	std::swap(query, match);

	queryNode = query.substr(0, query.size()-2); // remove leftmark and rightmark
	rightKmer = queryNode.substr(queryNode.size() - k, k);                       //last k characters of the query node
	leftKmer = queryNode.substr(0, k);                                           //first k characters of the query node
	queryLeftMark = query.substr(query.size() - 2, 1);
	queryRightMark = query.substr(query.size() - 1, 1);

	matchNode = match.substr(0, match.size()-2);
	matchLeftKmer = matchNode.substr(0, k);
	matchRightKmer = matchNode.substr(matchNode.size() - k, k);
	matchLeftMark = match.substr(match.size() - 2, 1);
	matchRightMark = match.substr(match.size() - 1, 1);

}

if (matchLeftMark != "1")
continue; // nothing to do here

if (rightKmer.compare(matchLeftKmer) == 0)
{
	string node = queryNode + matchNode.substr(k, matchNode.size() - k);
	it_gS->second.first = "-1";     // these are done being glued
	it_gS->second.second = "-1";

	remove(query); remove(match); // scrub the table

	if (queryLeftMark == "0" && matchRightMark == "0")
	{
		if (debug)
			cout << "Case 1 match and output: " << matchNode << " " << matchLeftMark << matchRightMark << "\n";
		output(node);
	}
	else    // at least one of left and right needs to be updated
	{
		insert(node, queryLeftMark == "1", matchRightMark == "1");
	}
}
else
{
	cout << "uh oh, not matching" << endl;
	cout << "  query: " << debug_highlight(query,rightKmer) <<  "\n";
	cout << "  match: " << debug_highlight(match,matchLeftKmer) << std::endl;
}
}

for (auto it = glueStrings.begin(); it != glueStrings.end();)
{
	if ((it->second).first == "-1" && (it->second).second == "-1")
	{
#ifdef SPARSEHASH
		glueStrings.erase(it);
#else
		it = glueStrings.erase(it);
#endif
	}
	else
		it ++;
}
#ifdef SPARSEHASH
glueStrings.resize(0); // effectively remove erased entries from memory
#endif
//cout << "new size: " << glueStrings.size() <<endl;



//	   cout << "what remains in the glue\n";
//	   for (auto it = glueStrings.begin(); it != glueStrings.end(); it++)
//	   {
//	   cout << debug_highlight((it->second).first,it->first) << " - " << debug_highlight((it->second).second,it->first) << "\n";
//	   {
//	   string seq = it->first;
//	   Model modelK1(kmerSize-1, minSize);
//	   Model::Kmer kmmerBegin=modelK1.codeSeed(seq.substr(0,kmerSize-1).c_str(),Data::ASCII);
//	   size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
//	   Model::Kmer kmmerEnd=modelK1.codeSeed(seq.substr(seq.size()-kmerSize+1,kmerSize-1).c_str(),Data::ASCII);
//	   size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));
//	   }
//	   }
//	  
}

*/

