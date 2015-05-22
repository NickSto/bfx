# This script takes a file with sequences (not fasta) separated by newlines and 
# counts the number of single, double, triple, quadruple and more than
# quintuple repeats of a particular sequence. It counts each repeat
# as one unit. In this case is the telomere hexamer sequence, but it can be 
# replaced by any sequence. It prints five files, one for each repeat count. August 2013 

#Author: Samarth Rangavittal, Rebeca Campos

# Usage: python findmotifseparate.py infile_name motif_search


import sys

infile = sys.argv[1]
fp_in = open(infile,'r')
motif = sys.argv[2]

single = open("single",'w')
double = open("double",'w')
triple = open("triple",'w')
quad = open("quad",'w')
above5 = open("above5",'w')

for line in fp_in :
	consider = str(line[:])
	iteration = len(consider) - len(motif) 
	single_count = 0
	double_count = 0
	triple_count = 0
	quad_count = 0
	above5_count = 0
	repeat_count = 0
	i = 0
	
	while i < iteration :
		suspect = consider[i:i+len(motif)]
		if suspect == motif :	
			while suspect == motif :
				i = i+len(motif)
				suspect = consider[i:i+len(motif)]
				repeat_count += 1
				
			if repeat_count == 1:
				single_count = single_count + 1
				repeat_count = 0
				
			elif repeat_count == 2:
				double_count = double_count + 1
				repeat_count = 0
				
			elif repeat_count == 3:
				triple_count = triple_count + 1
				repeat_count = 0
				
			elif repeat_count == 4:
				quad_count = quad_count + 1
				repeat_count = 0
			
			elif repeat_count >= 5:
				above5_count = above5_count + 1
				repeat_count = 0
		else :	
			i = i+1
			
	single.write(str(single_count)+"\n")
	double.write(str(double_count)+"\n")
	triple.write(str(triple_count)+"\n")
	quad.write(str(quad_count)+"\n")
	above5.write(str(above5_count)+"\n")
	
	#print "Single: ", single_count
	#print "Double: ", double_count
	#print "Triple: ", triple_count
	#print "Quad: ", quad_count
	#print "Above: ", above_count
	
fp_in.close()


