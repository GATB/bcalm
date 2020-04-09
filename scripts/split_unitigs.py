#!/usr/bin/env python
import sys, os
if len(sys.argv) < 4:
    print("split BCALM unitigs at reference extremities. More specifically:")
    print("  the script considers the sets B and E of all k-mers that are extremities of the reference genomes/contigs, respectively B for beginnings of contigs and E for ends.")
    print("  output is: modified unitigs such that each k-mer in B should be the beginning of an unitig, and each kmer in E should be end of an unitig.")
    print("  in order words, unitigs are split at kmers that are extremities of the reference sequences")
    print("This script is a small modification of pufferize.py")
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
# NOTE: unitigs links are discarded
output = open(unitigs+".split.fa","w")
nb_unitigs=-1
def save_unitig(header,seq):
    global output, p, nb_unitigs
    nb_unitigs += 1
    output.write(">unitig%s\n%s\n" % (nb_unitigs, seq))
    return(nb_unitigs)


unitig_skmer = {}
unitig_ekmer = {}

def create_unitig(header, unitig):
    global unitig_skmer, unitig_ekmer
    if len(unitig) == k:
        unitig = normalize(unitig)
    unitig_id=save_unitig(header, unitig)
    if normalize(unitig[:k]) in unitig_skmer:
        print("Warning, start kmer (%s) was also seen at start of unitig %s." % (unitig[:k],str(unitig_skmer[normalize(unitig[:k])])))
    if normalize(unitig[:k]) in unitig_ekmer:
        print("Warning, start kmer (%s) was also seen at end of unitig %s." % (unitig[:k],str(unitig_ekmer[normalize(unitig[:k])])))
    if normalize(unitig[-k:]) in unitig_skmer:
        print("Warning, last kmer (%s) was also seen at start of unitig %s." % (unitig[-k:],str(unitig_skmer[normalize(unitig[-k:])])))
    if normalize(unitig[-k:]) in unitig_ekmer:
        print("Warning, last kmer (%s) was also seen at end of unitig %s." % (unitig[-k:],str(unitig_ekmer[normalize(unitig[-k:])])))
    unitig_ekmer[normalize(unitig[-k:])] = [unitig_id, len(unitig)]
    unitig_skmer[normalize(unitig[:k])] = [unitig_id, len(unitig)]
   
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


output.close()
print("done. result is in: %s.split.fa" % unitigs)
