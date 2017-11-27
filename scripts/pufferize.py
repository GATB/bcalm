#!/usr/bin/env python
import sys, os
if len(sys.argv) < 4:
    print("alters BCALM's unitigs so that they fit pufferfish input. more specifically:")
    print("  the script considers the set B and E of all k-mers that are extremities of the reference genomes, respectively beginning and end of genomes.")
    print("  output is: modified unitigs such that each k-mer in B should be the beginning of an unitig, and each kmer in E should be end of an unitig.")
    print("  in order words, unitigs are split at kmers that are extremities of the reference sequences")
    exit("arguments: references.fa unitigs.fa k")

references=sys.argv[1]
unitigs=sys.argv[2]
k=int(sys.argv[3])

# https://www.biostars.org/p/710/#1412
from itertools import groupby
def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def revcomp(seq):
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def normalize(kmer):
    return kmer if kmer < revcomp(kmer) else revcomp(kmer)


#ref_kmers=set()
# parse references
#for header, ref in fasta_iter(references):
#  for kmer in [ref[:k], ref[-k:]]:
#        ref_kmers.add(normalize(kmer))
#        print(kmer)

# go through all reference strings and keep the starting and end kmers for each of them in the ref_skmers and ref_ekmers respectively.
ref_skmers=set()
ref_ekmers=set()
# parse references
for header, ref in fasta_iter(references):
    ref_skmers.add(ref[:k])
    ref_ekmers.add(ref[-k:])


# parse unitigs and split if necessary
# ASSUMPTION: we might need to split a unitig multiple times
# NOTE: the exact position of the split in unitig depends on whether the seen kmer is the first or last kmer in a reference string.
# NOTE: unitigs are renumbered consecutively
# NOTE: former unitigs links are discarded
#output = open(unitigs+".pufferized.fa","w")
output = open(unitigs+".pufferized.gfa", "w")
nb_unitigs=-1
def save_unitig(header,seq):
    global output, p, nb_unitigs
    nb_unitigs += 1
    output.write("S\t%s\t%s\n" % (nb_unitigs, seq))
    return(nb_unitigs)


unitig_skmer = {}
unitig_ekmer = {}

def create_unitig(header, unitig):
    global unitig_skmer, unitig_ekmer
    if len(unitig) == k:
        unitig = normalize(unitig)
    unitig_id=save_unitig(header, unitig)
    if normalize(unitig[:k]) in unitig_skmer or normalize(unitig[:k]) in unitig_ekmer:
        exit("Error: Initial kmer is repeated.")
    if normalize(unitig[-k:]) in unitig_skmer or normalize(unitig[-k:]) in unitig_ekmer:
        exit("Error: Last kmer is repeated.")
    unitig_ekmer[normalize(unitig[-k:])] = [unitig_id, len(unitig)]
   
print("Start parsing and spliting unitigs .. ")
for header, unitig in fasta_iter(unitigs):
    prev = 0
    for i in range(0,len(unitig)-k+1):
        kmer = unitig[i:i+k]       
        # cut up until first kmer but not the kmer itself
        if kmer in ref_skmers or revcomp(kmer) in ref_ekmers:
            if i+k-1-prev >= k:
                create_unitig(header, unitig[prev:i+k-1])
                prev = i
        # cut the unitig until the kmer, including it
        if kmer in ref_ekmers or revcomp(kmer) in ref_skmers:
            create_unitig(header, unitig[prev:i+k])
            prev = i+1
    #add the last and right most unitig:
    if len(unitig)-prev >= k:
        create_unitig(header, unitig[prev:])


print("Start reconstructing the path .. ")
pathCtr = 0
for header, ref in fasta_iter(references):
    output.write("\nP\t")
    i = 0   
    seq = ''
    while (i < len(ref)-k+1):
        kmer = ref[i:i+k]
        normalizedkmer = normalize(kmer)        
        unitigInfo = []
        ori = "+"
        if normalizedkmer in unitig_skmer and normalizedkmer in unitig_ekmer:
            unitigInfo = unitig_skmer[normalizedkmer]
            if kmer == normalizedkmer:
                ori = "+"
            else:
                ori = "-"
        elif normalizedkmer in unitig_skmer:
            unitigInfo = unitig_skmer[normalizedkmer]
            ori = "+"
        elif normalizedkmer in unitig_ekmer:
            unitigInfo = unitig_ekmer[normalizedkmer]
            ori = "-"
        else:
            print(pathCtr, " paths reconstructed.")
            exit("ERROR: kmer is not found in the start or end of a unitig \n{0}  ,  {1}".format(kmer, normalizedkmer))
        output.write("%s%s," % (unitigInfo[0], ori))
        #increase i by the total number of kmers in the current unitig
        i += (unitigInfo[1]-k+1)    
    pathCtr+=1

output.close()
print("done. result is in: %s.pufferized.gfa" % unitigs)
print("to update unitig links in the fasta header (necessary to get a GFA file), run:")
print("mv %s.pufferized.fa %s" % (unitigs,unitigs))
prefix = "[prefix]" if ("unitigs.fa" not in unitigs) else (unitigs.split(".unitigs.fa")[0]+".h5")
script_dir = os.path.dirname(os.path.realpath(__file__))
bcalm_path = "bcalm" if not os.path.isfile("%s/../build/bcalm" % script_dir) else os.path.abspath("%s/../build/bcalm" %script_dir)
print("%s -in %s -skip-bcalm -skip-bglue -redo-links" % (bcalm_path,prefix))
print("%s/convertToGFA.py %s %s.gfa %d" % (script_dir,unitigs,unitigs,k))
