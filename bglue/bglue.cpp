/* remaining issue: 
 * - inconsistency of results between machines (laptop lifl vs cyberstar231) */
#include <gatb/gatb_core.hpp>
#include "unionFind.hpp"
#include "glue.hpp"
#include <atomic>

class bglue : public Tool
{
public:
    bglue() : Tool ("bglue"){
	getParser()->push_back (new OptionOneParam ("-k", "kmer size",  false,"31"));
	getParser()->push_back (new OptionOneParam ("-in", "input file",  true)); // necessary for repartitor
	getParser()->push_back (new OptionOneParam ("-out", "output file",  false, "out.fa"));
	getParser()->push_front (new OptionNoParam  ("--only-uf",   "(for debugging only) stop after UF construction", false));
    };

    // Actual job done by the tool is here
    void execute ();
};

using namespace std;

const size_t SPAN = KMER_SPAN(1); // TODO: adapt span using Minia's technique
typedef Kmer<SPAN>::Type  Type;
typedef Kmer<SPAN>::Count Count;
typedef Kmer<SPAN>::ModelCanonical ModelCanon;
typedef Kmer<SPAN>::ModelMinimizer <ModelCanon> Model;
size_t kmerSize=31;
size_t minSize=8;

    // a hash wrapper for hashing kmers
    template <typename ModelType>
    class Hasher_T 
    {
       public:
        ModelType model;
        
        Hasher_T(ModelType &model) : model(model) {};

        uint64_t operator ()  (const typename ModelType::Kmer& key, uint64_t seed = 0) const  {
                return model.getHash(key.value());
                }
    };


extern char revcomp (char s); // glue.cpp

string rc(string &s) {
	string rc;
	for (signed int i = s.length() - 1; i >= 0; i--) {rc += revcomp(s[i]);}
	return rc;
}


struct markedSeq
{
    string seq;
    bool lmark, rmark;
    string ks, ke; // [start,end] kmers of seq, in canonical form (redundant information with seq, but helpful)

    markedSeq(string seq, bool lmark, bool rmark, string ks, string ke) : seq(seq), lmark(lmark), rmark(rmark), ks(ks), ke(ke) {};

    void revcomp()
    {
        seq = rc(seq);
        std::swap(lmark, rmark);
        std::swap(ks, ke);
    }
};

vector<vector<markedSeq> > determine_order_sequences(vector<markedSeq> sequences)
{
    bool debug = false ;
    unordered_map<string, set<int> > kmerIndex;
    set<int> usedSeq;
    vector<vector<markedSeq>> res;
    unsigned int nb_chained = 0;

    // index kmers to their seq
    for (unsigned int i = 0; i < sequences.size(); i++)
    {
        kmerIndex[sequences[i].ks].insert(i);
        kmerIndex[sequences[i].ke].insert(i);
    }

    for (unsigned int i = 0; i < sequences.size(); i++)
    {
        markedSeq current = sequences[i];
        if (usedSeq.find(i) != usedSeq.end())
                continue; // this sequence has already been glued

        if (current.lmark & current.rmark)
            continue; // not the extremity of a chain

        if (current.lmark)
            current.revcomp(); // reverse so that lmark is false

        assert(current.lmark == false);

        vector<markedSeq> chain;
        chain.push_back(current);
    
        bool rmark = current.rmark;
        int current_index = i;
        int start_index = i;
        usedSeq.insert(i);

        while (rmark)
        {
            if (debug)
                std::cout << "current ke " << current.ke << " index " << current_index << " markings: " << current.lmark << current.rmark <<std::endl;

            set<int> candidateSuccessors = kmerIndex[current.ke];
            assert(candidateSuccessors.find(current_index) != candidateSuccessors.end());
            candidateSuccessors.erase(current_index);
            assert(candidateSuccessors.size() == 1); // normally there is exactly one sequence to glue with

            int successor_index = *candidateSuccessors.begin(); // pop()
            assert(successor_index != current_index);
            markedSeq successor = sequences[successor_index]; 

            if (successor.ks != current.ke || (!successor.lmark))
                successor.revcomp();
            
            if (debug)
                std::cout << "successor " << successor_index << " successor ks ke "  << successor.ks << " "<< successor.ke << " markings: " << successor.lmark << successor.rmark << std::endl;
            
            assert(successor.lmark);
            assert(successor.ks == current.ke);

            if (successor.ks == successor.ke)
            {
                if (debug)
                    std::cout << "successor seq loops: " << successor.seq << std::endl;
                assert(successor.seq.size() == kmerSize);
                if (successor.lmark == false)
                    assert(successor.rmark == true); 
                else
                    assert(successor.rmark == false);
                // it's the only possible cases I can think of
                
                // there is actually nothing to be done now, it's an extremity, so it will end.
                // on a side note, it's pointless to save this kmer in bcalm.
            }
            

            current = successor;
            chain.push_back(current);
            current_index = successor_index;
            rmark = current.rmark;
            assert((usedSeq.find(current_index) == usedSeq.end()));
            usedSeq.insert(current_index);

        }
        res.push_back(chain);
        nb_chained += chain.size();
    }
    assert(sequences.size() == nb_chained); // make sure we've scheduled to glue all sequences in this partition
    return res;
}

/* straightforward glueing of a chain
 * sequences should be ordered and in the right orientation
 * so, it' just a matter of chopping of the first kmer
 */
string glue_sequences(vector<markedSeq> chain)
{
    string res;
    string previous_kmer = "";
    unsigned int k = kmerSize;
    bool last_rmark = false;

    for (auto it = chain.begin(); it != chain.end(); it++)
    {
        string seq = it->seq;

        if (previous_kmer.size() == 0) // it's the first element in a chain
        {
            assert(it->lmark == false);
            res += seq;
        }
        else
        {
            assert(seq.substr(0, k).compare(previous_kmer) == 0);
            res += seq.substr(k);
        }
        
        previous_kmer = seq.substr(seq.size() - k);
        assert(previous_kmer.size() == k);
        last_rmark = it->rmark;
    }
    assert(last_rmark == false);
    if (last_rmark) { cout<<"bad gluing, missed an element" << endl; exit(1); } // in case assert()'s are disabled

    return res;
}

void output(string seq, BankFasta &out)
{
    Sequence s (Data::ASCII);
    s.getData().setRef ((char*)seq.c_str(), seq.size());
    out.insert(s);
}


void bglue::execute (){

    kmerSize=getInput()->getInt("-k");
    size_t k = kmerSize;
    string inputFile(getInput()->getStr("-in")); // necessary for repartitor

    string h5_prefix = inputFile.substr(0,inputFile.size()-2);
    BankFasta in (h5_prefix + "glue");


    Storage* storage = StorageFactory(STORAGE_HDF5).load ( inputFile.c_str() );

    LOCAL (storage);
    /** We get the dsk and minimizers hash group in the storage object. */
    Group& dskGroup = storage->getGroup("dsk");
    Group& minimizersGroup = storage->getGroup("minimizers");

    typedef Kmer<SPAN>::Count Count;
    Partition<Count>& partition = dskGroup.getPartition<Count> ("solid");
    size_t nbPartitions = partition.size();
    cout << "DSK created " << nbPartitions << " partitions" << endl;

    /** We retrieve the minimizers distribution from the solid kmers storage. */
    Repartitor repart;
    repart.load (minimizersGroup);

    u_int64_t rg = ((u_int64_t)1 << (2*minSize));

    /* Retrieve frequency of minimizers;
     * actually only used in minimizerMin and minimizerMax */
    uint32_t *freq_order = NULL;

    int minimizer_type = 1; // typical bcalm usage.
    if (minimizer_type == 1) 
    {
        freq_order = new uint32_t[rg];
        Storage::istream is (minimizersGroup, "minimFrequency");
        is.read ((char*)freq_order, sizeof(uint32_t) * rg);
    }

    
    Model model(kmerSize, minSize, Kmer<SPAN>::ComparatorMinimizerFrequencyOrLex(), freq_order);
    Model modelK1(kmerSize-1, minSize,  Kmer<SPAN>::ComparatorMinimizerFrequencyOrLex(), freq_order);
    Model modelK2(kmerSize-2, minSize,  Kmer<SPAN>::ComparatorMinimizerFrequencyOrLex(), freq_order);
    Model modelM(minSize, minSize,  Kmer<SPAN>::ComparatorMinimizerFrequencyOrLex(), freq_order);

    // create a hasher for UF
    ModelCanon modelCanon(kmerSize); // i'm a bit lost with those models.. I think GATB could be made more simple here.
    Hasher_T<ModelCanon> hasher(modelCanon);

    // create a UF data structure
#if 0
    unionFind<unsigned int> ufmin;
    unionFind<std::string> ufprefixes; 
    unsigned int prefix_length = 10;
    unionFind<std::string> ufkmerstr; 
#endif
    unionFind<std::string> ufkmerstr(100); 
    // those were toy one, here is the real one:

    typedef uint32_t partition_t;
    unionFind<partition_t> ufkmers(1000); 
    // instead of UF of kmers, do a union find of hashes of kmers. less memory. will have collisions, but that's okay i think. let's see.
    // also, use many locks

    
    // We create an iterator over this bank.
    BankFasta::Iterator it (in);

    // We loop over sequences.
    /*for (it.first(); !it.isDone(); it.next())
    {
        string seq = it->toString();*/
    auto createUF = [k, &modelK1 /* crashes if copied!*/, &model, \
        &ufkmerstr, &ufkmers, &hasher](const Sequence& sequence)
    {
        string seq = sequence.toString();

        if (seq.size() < k) 
        {
            std::cout << "unexpectedly small sequence found ("<<seq.size()<<"). did you set k correctly?" <<std::endl; exit(1);
        }
        string comment = sequence.getComment();
        bool lmark = comment[0] == '1';
        bool rmark = comment[1] == '1';

        if (!(lmark & rmark)) // if either mark is 0, it's an unitig extremity, so nothing to glue here
            return;

        string kmerBegin = seq.substr(0, k );
        string kmerEnd = seq.substr(seq.size() - k , k );

        Model::Kmer kmmerBegin = model.codeSeed(kmerBegin.c_str(), Data::ASCII);
        Model::Kmer kmmerEnd = model.codeSeed(kmerEnd.c_str(), Data::ASCII);
      
        ufkmers.union_(hasher(kmmerBegin), hasher(kmmerEnd));

#if 0 
        string canonicalKmerBegin = modelK1.toString(kmmerBegin.value());
        string canonicalKmerEnd = modelK1.toString(kmmerEnd.value());
        ufkmerstr.union_(canonicalKmerBegin, canonicalKmerEnd);

        size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
        size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));
        ufmin.union_(leftMin, rightMin);

        string prefixCanonicalKmerBegin = canonicalKmerBegin.substr(0, prefix_length);
        string prefixCanonicalKmerEnd = canonicalKmerEnd.substr(0, prefix_length);
        ufprefixes.union_(prefixCanonicalKmerBegin, prefixCanonicalKmerEnd);
#endif    
        

    };

    //setDispatcher (new SerialDispatcher()); // force single thread
    setDispatcher (  new Dispatcher (getInput()->getInt(STR_NB_CORES)) );
    getDispatcher()->iterate (it, createUF);

#if 0
    ufmin.printStats("uf minimizers");

    ufprefixes.printStats("uf " + to_string(prefix_length) + "-prefixes of kmers");
    
    ufkmerstr.printStats("uf kmers, std::string");
#endif
    

    if (getParser()->saw("--only-uf")) // for debugging
        return;

    ufkmers.printStats("uf kmers");

    // now glue using the UF
    BankFasta out (getInput()->getStr("-out"));

    auto decide_if_process_seq_in_partition = [](bool found_partition, uint32_t partition, uint32_t pass, uint32_t nb_passes)
    {
        if (!found_partition)
        {   
            return (pass == (nb_passes - 1)); // no-partition sequences (those who do need to be glued) are output at the very end.
        }

        return (partition % nb_passes) == pass;
    };

    int nb_passes = 4; // TODO: tune that

    // to reduce memory usage, process the glued file in multiple passes.
    for (int pass = 0; pass < nb_passes; pass++)
    {
        cout << "pass " << pass << endl;
        unordered_map<int,vector<markedSeq>> msInPart;
        for (it.first(); !it.isDone(); it.next()) // I don't think this loop could be parallelized at all. if it did, all threads would compute the same partition index.
            // par contre; pourrait y avoir une iteration en parallele (avec un autre set de threads) qui calculeraient l'index de partition.
        {
            string seq = it->toString();

            string kmerBegin = seq.substr(0, k );
            string kmerEnd = seq.substr(seq.size() - k , k );

            // make canonical kmer
            Model::Kmer kmmerBegin = model.codeSeed(kmerBegin.c_str(), Data::ASCII);
            Model::Kmer kmmerEnd = model.codeSeed(kmerEnd.c_str(), Data::ASCII);

            partition_t partition = 0;
            bool found_partition = false;

            string comment = it->getComment();
            bool lmark = comment[0] == '1';
            bool rmark = comment[1] == '1';


            if (lmark)
            {
                found_partition = true;
                partition = ufkmers.find(hasher(kmmerBegin));
            }

            if (rmark)
            {
                if (found_partition) // just do a small check
                {   
                    if (ufkmers.find(hasher(kmmerEnd)) != partition)
                    { std::cout << "bad UF! left kmer has partition " << partition << " but right kmer has partition " << ufkmers.find(hasher(kmmerEnd)) << std::endl; exit(1); }
                }

                partition = ufkmers.find(hasher(kmmerEnd));
                found_partition = true;
            }

            // at this point, decide if we process this sequence in this pass or not
            if (!decide_if_process_seq_in_partition(found_partition, partition, pass, nb_passes))
                continue;

            string ks = model.toString(kmmerBegin.value());
            string ke = model.toString(kmmerEnd  .value());
            markedSeq ms(seq, lmark, rmark, ks, ke);

            if (!found_partition) // this one doesn't need to be glued 
            {
                output(ms.seq, out);
                continue;
            }

            msInPart[partition].push_back(ms);
        }


        // now iterates all sequences in a partition to glue them in clever order (avoid intermediate gluing)
        for (auto it = msInPart.begin(); it != msInPart.end(); it++)
        {
            // TODO: do a task based parallelism here!!! perfect opportunity.

            //std::cout << "1.processing partition " << it->first << std::endl;
            vector<vector<markedSeq>> ordered_sequences = determine_order_sequences(it->second);
            //std::cout << "2.processing partition " << it->first << " nb ordered sequences: " << ordered_sequences.size() << std::endl;

            for (auto itO = ordered_sequences.begin(); itO != ordered_sequences.end(); itO++)
            {
                string seq = glue_sequences(*itO);

                // TODO: use a output queue, or a lock, or a task based output again, i don't know
                output(seq, out);
            }
        }
    }

//#define ORIGINAL_GLUE
#ifdef ORIGINAL_GLUE
    // We loop again over sequences
    // but this time we glue!
    int nb_glues = 1;
    BankFasta out (getInput()->getStr("-out") + ".original_glue_version");
    GlueCommander glue_commander(kmerSize, &out, nb_glues, &model);
    for (it.first(); !it.isDone(); it.next())
    {
        string seq = it->toString();
        string comment = it->getComment();
        bool lmark = comment[0] == '1';
        bool rmark = comment[1] == '1';
        GlueEntry e(seq, lmark, rmark, kmerSize);
        glue_commander.insert(e);
    }

    glue_commander.stop();
    cout << "Final glue:\n";
    glue_commander.dump();
    cout << "*****\n";
    glue_commander.printMemStats();
#endif

}

int main (int argc, char* argv[])
{
    try
    {
        bglue().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

