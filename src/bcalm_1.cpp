#include <gatb/debruijn/impl/GraphUnitigs.hpp>
#include <bcalm_1.hpp>

using namespace std;

/*
 * where did all the code go? now it's mostly in ../gatb-core/gatb-core/src/gatb/bcalm/
 */

/********************************************************************************/


bcalm_1::bcalm_1 ()  : Tool ("bcalm_1"){
    // old options, now using GATB's built-in options for kmer counting and graph creation
    // but TODO would be nice to integrate --nb-glue-partitions in gatb someday
/*	getParser()->push_back (new OptionOneParam ("-in", "input file",  true));
	getParser()->push_back (new OptionOneParam ("-out", "output prefix",  false, "unitigs"));
	getParser()->push_back (new OptionOneParam ("-k", "kmer size",  false,"31"));
	getParser()->push_back (new OptionOneParam ("-m", "minimizer size",  false,"8"));
	getParser()->push_back (new OptionOneParam ("-abundance", "abundance threshold",  false,"1"));
	getParser()->push_back (new OptionOneParam ("-minimizer-type", "use lexicographical minimizers (0) or frequency based (1)",  false,"1"));
	getParser()->push_back (new OptionOneParam ("-dsk-memory", "max memory for kmer counting (MB)", false, "1500"));
	getParser()->push_back (new OptionOneParam ("-dsk-disk", "max disk space for kmer counting (MB)", false, "default"));

    // glue options
    getParser()->push_front (new OptionNoParam  ("--only-uf",   "(for debugging only) stop after UF construction", false));
	getParser()->push_front (new OptionNoParam  ("--uf-stats",   "display UF statistics", false));
	getParser()->push_back (new OptionOneParam ("--nb-glue-partitions", "number of glue files on disk",  false,"200"));
*/

    IOptionsParser* graphParser = GraphUnitigsTemplate<32>::getOptionsParser(false);

    // hiding options
    if (IOptionsParser* p = graphParser->getParser(STR_KMER_ABUNDANCE_MIN_THRESHOLD))  {  p->setVisible(false); }
    if (IOptionsParser* p = graphParser->getParser(STR_HISTOGRAM_MAX))  {  p->setVisible(false); }
    if (IOptionsParser* p = graphParser->getParser(STR_SOLIDITY_KIND))  {  p->setVisible(false); } // oohh. multi-sample dbg construction someday maybe?
    if (IOptionsParser* p = graphParser->getParser(STR_URI_SOLID_KMERS))  {  p->setVisible(false); }
    
    // setting defaults
    if (Option* p = dynamic_cast<Option*> (graphParser->getParser(STR_REPARTITION_TYPE)))  {  p->setDefaultValue ("1"); }
    if (Option* p = dynamic_cast<Option*> (graphParser->getParser(STR_MINIMIZER_TYPE)))  {  p->setDefaultValue ("1"); }

    getParser()->push_back(graphParser);


}


template <size_t span>
struct Functor  {  void operator ()  (bcalm_1 *bcalm)
    {
        typedef GraphUnitigsTemplate<span> GraphType;
        GraphType graph;

        if (bcalm->getInput()->get(STR_URI_INPUT) != 0)
        {
            graph = GraphType::create (bcalm->getInput(), false /* do not load unitigs after*/);
        }
        else
        {
            throw OptionFailure (bcalm->getParser(), "Specifiy -in");
        }

        // delete the .h5 file
        bool delete_h5_file = true;
        if (delete_h5_file)
        {
            // copies the h5 naming mechanism in GraphUnitigs.cpp
            string input = bcalm->getInput()->getStr(STR_URI_INPUT);
            string prefix;
            if (bcalm->getInput()->get(STR_URI_OUTPUT)) 
                prefix = bcalm->getInput()->getStr(STR_URI_OUTPUT);
            else
                prefix = System::file().getBaseName (input) + prefix;

            System::file().remove (prefix + ".h5");
        }

    }
};

void bcalm_1::execute (){

    std::cout << "BCALM 2, version " << VERSION; 
#ifdef GIT_SHA1
    std::cout << ", git commit " << GIT_SHA1;
#endif
    std::cout << std::endl;


    /** we get the kmer size chosen by the end user. */
    size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

    /** We launch Minia with the correct Integer implementation according to the choosen kmer size. */
    Integer::apply<Functor,bcalm_1*> (kmerSize, this);

}


