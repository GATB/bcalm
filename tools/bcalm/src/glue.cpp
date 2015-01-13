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

/*
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
*/

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

string rcnorm ( string seq ) {
	return std::min (seq, reversecomplement(seq));
}

char bits2char (bool leftbit, bool rightbit) {
	if (leftbit == 0 && rightbit == 0) return 'A';
	else if (leftbit == 0 && rightbit == 1) return 'C';
	else if (leftbit == 1 && rightbit == 0) return 'G';
	else return 'T';
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


string tostring(GlueEntry e, string key) {
	ostringstream out;
	if (e.getLmark()) out << "_"; else out << " ";
	out << debug_highlight(e.getSeq(), key);
	if (e.getRmark()) out << "_"; else out << " ";
	return out.str();
}

GlueEntry::GlueEntry(RawEntry raw, size_t _kmerSize){
	kmerSize = _kmerSize;
	seq = "";
#ifdef TWOBITGLUEHASH

	//convert encoded string raw to a binary string bits
	RawEntry bits;
	for (int i = 0; i < raw.size(); i++) {
		bitset<8> byte(int(raw[i]));
		string bits_aux = byte.to_string();
		transform(bits_aux.begin(), bits_aux.end(), bits_aux.begin(),
				          bind2nd(std::minus<int>(), '0')); 
		bits.insert(bits.end(), bits_aux.begin(), bits_aux.end()); // += bits_aux2;
	}
	//cout << "Constructor extracted bits: \nbits: " << bits << ".\n";
	if (bits.size() != raw.size() * 8) {
		cout << "fail in GlueEntry::GlueEntry bit conversion.\n";
		exit(1);
	}

	lmark = bits[0];
	rmark = bits[1];
	
	seq = "";
	int leftoverbits = 2 * ( 2 * bits[2] + bits[3]); //to represent how many bits are to be read from the last byte: 00 = 0 bits, 01 = 2 bits, 10 = 4 bits, 11 = 6 bits
	if (leftoverbits == 0) leftoverbits = 8;
	for (int i = 4; i < bits.size() - 8 + leftoverbits; i += 2) {
		seq.push_back(bits2char(bits[i], bits[i+1]));
	}
#else
	seq = raw.substr(0, raw.size() - 2);
	lmark = raw.at(raw.length() - 2);
	rmark = raw.at(raw.length() - 1);
#endif
}


uint8_t bits2byte (RawEntry bits, int start) {
	return  
	1 * bits[start + 7] +  
	2 * bits[start + 6] + 
	4 * bits[start + 5] + 
	8 * bits[start + 4] + 
	16 * bits[start + 3] + 
	32 * bits[start + 2] + 
	64 * bits[start + 1] + 
	128 * bits[start];
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
	RawEntry bits(rawLengthBytes * 8, 0);
	bits[0] = int(lmark);
	bits[1] = int(rmark);
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
	//string raw(rawLengthBytes, 0);
	RawEntry raw;
	for (int i = 0; i < bits.size(); i += 8) {
		raw.push_back(bits2byte(bits, i));
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
	
	if (!match.getLmark()  || (query.getRkmer() != match.getLkmer()) || (query.getRkmer() != key)) {
		std::swap(query, match);
	}

	if (!(match.getLmark() && (query.getRkmer() == match.getLkmer()) && (query.getRkmer() == key)) ) {
		cout << "glueSingleEntry received non-gluable input:\n";
		cout << "  query: " << tostring(query, key) << endl;
		cout << "  match: " << tostring(match, key) << endl;
		exit(1);
	}

	string gluedStr = query.getSeq() + match.getSeq().substr(kmerSize, match.getSeq().size() - kmerSize);
	glueResult = GlueEntry(gluedStr, query.getLmark(), match.getRmark(), kmerSize);
}


bool Glue::check_if_empty(GlueEntry newEntry, string key) {
	GlueEntry e;

	if (key.compare(reversecomplement(key)) >= 0) {
		newEntry.revComp();
		key = reversecomplement(key);
	}

	if (!glueStorage.find(key, e)) {
		return true;
	} else if (e.isEmpty()) {
		return true;
	}
	return false;
}


/*void Glue::write2hash(string key, GlueEntry e) {
	//auto tag = tagSystem.createTag(e);
	//glueStorage.insertAtKey(key, tag);
	glueStorage.insertAtKey(key, e);
}
*/

void Glue::insert_aux(GlueEntry newEntry, string key, GlueEntry & glueResult) { 
	GlueEntry e;

	if (key.compare(reversecomplement(key)) >= 0) {
		newEntry.revComp();
		key = reversecomplement(key);
	}

	//if (key == key_of_interest)  cout << key << "\tinsert_aux\t" << tostring(newEntry, key)  << "\tprior\t" << glueStorage.dump(key, false) << "\t";
	if (glueDebug) cout << key << "\tinsert_aux\t" << tostring(newEntry, key)  << "\tprior\t" << glueStorage.dump(key, false) << "\t";

	if (!glueStorage.find(key, e)) {
		glueStorage.insertAtKey(key, newEntry);
	} else {
		if (e.isEmpty()) {
			glueStorage.insertAtKey(key, newEntry);
		} else {
			//Looks like we need to do a glue!
			glueResult = GlueEntry();
			glueSingleEntry(e, newEntry, key, glueResult); //glue the two strings
			glueStorage.insertAtKey(key, GlueEntry()); //clear entry
		}
	}
	//if (key == key_of_interest) cout << "after\t" << glueStorage.dump(key,false ) << endl;
	if (glueDebug) cout << "after\t" << glueStorage.dump(key,false ) << endl;

}


void Glue::insert(GlueEntry e, bool process) {
	
	/* testing code 
	string str;
	while (cin >> str) {
		GlueEntry en(str, false, false, 2);
		RawEntry raw = en.getRaw();
		cout << tostring(en, "Z") << endl << raw << endl << tostring(GlueEntry(raw,2), "Z") << endl;
	}
	*/

	startTimer();
	bool oldGlueDebug = glueDebug;
	if (process == false) {
		cout << "insertion without processing is no longer supported\n";
		exit(1);
	}

	
	if ((e.getSeq().find(key_of_interest) != std::string::npos) || (e.getSeq().find(reversecomplement(key_of_interest)) != std::string::npos)) {
		glueDebug = true;
	}

	if (glueDebug) 
		cout << "insert\t" << tostring(e,"") << endl;

	GlueEntry insRes;

	if (!e.getRmark() && !e.getLmark()) {
		output(e.getSeq());
	} else if (e.getLmark() && !e.getRmark()) {
		insert_aux(e, e.getLkmer(), insRes);
	} else if (!e.getLmark() && e.getRmark()) {
		insert_aux(e, e.getRkmer(), insRes);
	} else { 
		//pick one that is not empty
		//this is not necessary for correctness, picking an arbitrary one will do
		//but picking a non-empty one may force glues to happen sooner rather later, preventing build-ups of long sequential chains
		if (!check_if_empty(e, e.getRkmer())) {
			insert_aux(e, e.getRkmer(), insRes);
		} else {
			insert_aux(e, e.getLkmer(), insRes);
		}
	}
	

	if (!insRes.isEmpty()) insert(insRes, true);

	glueDebug = oldGlueDebug; 

	stopTimer();
}

void GlueStorage::cleanup() {
	//GlueEntry e("bla", false, true, kmerSize);
	GlueEntry e;
	for (auto it = glueMap.begin(); it != glueMap.end(); ) {
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

void Glue::output(string seq)
{
	Sequence s (Data::ASCII);
	s.getData().setRef ((char*)seq.c_str(), seq.size());
	out.insert(s);
}

void Glue::stopTimer() { 
	if (--timerReferenceCount == 0) {
		auto endTime = chrono::system_clock::now();
		totalTime += chrono::duration_cast<chrono::seconds>(endTime - startTime);
	}
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

