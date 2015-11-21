#include <gatb/gatb_core.hpp>
#include "unionFind.hpp"

class bglue : public Tool
{
public:
    bglue() : Tool ("bglue"){
	getParser()->push_back (new OptionOneParam ("-k", "kmer size",  false,"31"));
	getParser()->push_back (new OptionOneParam ("-in", "input file",  true)); // necessary for repartitor
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

void bglue::execute (){

    kmerSize=getInput()->getInt("-k");
    size_t k = kmerSize;
    string inputFile(getInput()->getStr("-in")); // necessary for repartitor

    BankFasta in ("out_to_glue");


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

    // create a UF data structure
    unionFind<unsigned int> ufmin;
    unionFind<std::string> ufkmers; 
    unionFind<std::string> ufprefixes; 
    unsigned int prefix_length = 10;
    
    // We create an iterator over this bank.
    BankFasta::Iterator it (in);

    unsigned long nbSmall = 0;
    // We loop over sequences.
    for (it.first(); !it.isDone(); it.next())
    {
        string seq = it->toString();

        if (seq.size() < k) 
        {
            nbSmall++;
            continue;
        }


        string kmerBegin = seq.substr(0, k - 1);
        string kmerEnd = seq.substr(seq.size() - k + 1, k - 1);

        Model::Kmer kmmerBegin = modelK1.codeSeed(kmerBegin.c_str(), Data::ASCII);
        size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
        Model::Kmer kmmerEnd = modelK1.codeSeed(kmerEnd.c_str(), Data::ASCII);
        size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));
       
        ufmin.union_(leftMin, rightMin);
       
        string canonicalKmerBegin = modelK1.toString(kmmerBegin.value());
        string canonicalKmerEnd = modelK1.toString(kmmerEnd.value());

        ufkmers.union_(canonicalKmerBegin, canonicalKmerEnd);
        
        string prefixCanonicalKmerBegin = canonicalKmerBegin.substr(0, prefix_length);
        string prefixCanonicalKmerEnd = canonicalKmerEnd.substr(0, prefix_length);

        ufprefixes.union_(prefixCanonicalKmerBegin, prefixCanonicalKmerEnd);
    }

    if (nbSmall)
        std::cout<< nbSmall << " compacted unitig were smaller than k?! " << std::endl;

    ufmin.printStats("uf minimizers");

    ufkmers.printStats("uf kmers");
   
    ufprefixes.printStats("uf " + to_string(prefix_length) + "-prefixes of kmers");
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

