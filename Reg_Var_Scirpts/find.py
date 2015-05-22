# February 20, 2014
# This script finds all possible 13-mer instability motifs (McVean) from Myers 2005
# and counts them (col1) in each line of DNA sequences (normally these will correspond to
# windows in the genome). 

# Author: Rebeca Campos

# Usage python find.py infile_name

import sys
import re


infile = sys.argv[1]
fp_in = open(infile,'r')

output= open('13-mer_instability','w')


for line in fp_in :
    consider = str(line[:])
    runs = re.findall(r"CC[ACGT]CC[ATGC]T[ATGC][ATCG]CC[ATGC]C", consider)
    #print runs
    length = len(runs)
    #print length
    motif_length = length * 13
    #print motif_length
    seq_length = len(consider)
    #print seq_length
    expected_13mer = float(seq_length)/13
    print expected_13mer
    ratio = float(motif_length)/seq_length
    #print ratio
    output.write(str(length)+"\t"+ str(motif_length)+"\t"+ str(seq_length)+"\t"+ str(ratio)+"\n")

fp_in.close()
