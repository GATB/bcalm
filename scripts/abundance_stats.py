
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

#!/usr/bin/env python
import sys, os
if len(sys.argv) < 2:
    print("prints some abundance statitics of a unitigs FASTA file produced by BCALM")
    exit("arguments: unitigs.fa")
 
unitigs = sys.argv[1]
abundances = []
<<<<<<< HEAD
from collections import defaultdict
totsize = defaultdict(int)
=======
>>>>>>> be8238dbac264fc60739c77effea906e3f01a8e5
for header, unitig in fasta_iter(unitigs):
    for field in header.split():
        if field.startswith("km:f:"):
            abundance = field.split(":")[-1]
            #print(abundance)
<<<<<<< HEAD
            abundance = int(float(abundance)) # convert to rounded int
            abundances += [abundance]
            totsize[abundance] += len(unitig)

from collections import Counter
c = Counter(abundances)
print("'value' : 'number of unitigs having this mean abundance value' : 'total size of unitigs having this mean abundance'")
for val in sorted(list(c)):
    print(val,":",c[val],':',totsize[val])
=======
            abundances += [int(float(abundance))]

from collections import Counter
c = Counter(abundances)
print("'value' : 'number of unitigs having this mean abundance value'")
for val in sorted(list(c)):
    print(val,":",c[val])
>>>>>>> be8238dbac264fc60739c77effea906e3f01a8e5
