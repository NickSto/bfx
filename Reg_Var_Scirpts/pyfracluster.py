#!/usr/bin/env python
"""
Seeks to find intervals for which the fraction of certain types of sites
exceeds a threshold.  The intervals found are such that *all* prefixes and
suffixes of an interval satisfy that threshold.

An interval with P positive and X neutral sites satifies the threshold F
if P/(P+X) >= F

		P/(P+X) >= F   =>   P >= FP+FX   =>   (1-F)P - FX >= 0

So we can identify such regions by keeping a running sum of (1-F)P - FX.
The latter is accomplished by adding 1-F to the sum whenever we encounter a
positive site, and subtracting F whenever we encounter a neutral site.

So we can use weigths of 1-F and -F for postive and negative sites,
respectively.

Note that it is possible to incorporate negative (or penalty) sites into this
paradigm, but we do not currently support those.

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys


def main():
	global debug

	rate      = None
	positives = "+"
	neutrals  = "-"
	blockers  = None
	minLength = 1
	origin    = "zero"
	debug     = []

	for arg in sys.argv[1:]:
		if (arg.startswith("--length=")):
			minLength = int(arg.split("=",1)[1])
		elif (arg.startswith("--origin=")):
			origin = arg.split("=",1)[1]
			if   (origin == "0"): origin = "zero"
			elif (origin == "1"): origin = "one"
			assert (origin in ["zero","one"])
		elif (arg.startswith("--blockers=")):
			blockers = arg.split("=",1)[1]
		elif (arg.startswith("{")) and (arg.endswith("}")):
			(positives,neutrals) = parse_spec(arg)		
		elif (arg.startswith("--debug=")):
			debug += (arg.split("=",1)[1]).split(",")
		elif (arg.startswith("--")):
			assert (False), "can't understand %s" % arg
		elif (rate == None):
			rate = float_or_fraction(arg)
		else:
			assert (False), "can't understand %s" % arg

	assert (rate != None), \
	       "you gotta give me a rate fraction"

	# process the sequences

	lineNumber = 0
	for line in sys.stdin:
		lineNumber += 1
		line = line.strip()
		fields = line.split()

		if (len(fields) == 1):
			name = "line_%d" % lineNumber
			seq  = line
		else:
			name = fields[0]
			seq  = " ".join(fields[1:])

		for (start,end) in dense_intervals(seq,rate,positives,neutrals,
		                                   blockers=blockers,
		                                   minLength=minLength):
			numer = 0
			for ix in range(start,end):
				if (seq[ix] in positives): numer += 1
			denom = end - start
			if (origin == "one"): start += 1
			print "%s\t%s\t%s\t%d/%d\t%s" \
			    % (name,start,end,numer,denom,seq[start:end])


# parse a string of the form {positives}/{positives_and_neutrals}

def parse_spec(s):
	if ("/" not in s): raise ValueError
	(n,d) = s.split("/",1)
	if (not n.startswith("{")) or (not n.endswith("}")): raise ValueError
	if (not d.startswith("{")) or (not d.endswith("}")): raise ValueError

	positives = n[1:-1]
	d         = d[1:-1]

	for ch in positives:
		if (ch not in d): raise ValueError

	neutrals = [ch for ch in d if (ch not in positives)]
	return (positives,neutrals)


# convert a string to a number, allowing fractions

def float_or_fraction(s):
	if ("/" in s):
		(numer,denom) = s.split("/",1)
		return float(numer)/float(denom)
	else:
		return float(s)


# dense_intervals--
#	Find all non-overlapping runs with a good enough rate (of positives), and
#	which meet our length threshold.
#
#	The algorithm used is adapted from Zhang, Berman, Miller, "Post-processing
#	long pairwise alignments", Bioinformatics Vol. 15 no. 12 1999.
#
# $$$ we use the denominator as the threshold, but we really should use the
# $$$ .. numerator, comparing it to minLength*rate

def dense_intervals(seq,rate,positives,neutrals,blockers="",minLength=None):

	if (blockers == None):
		blockers = "".join([chr(n) for n in range(1,256)
		                           if  (chr(n) not in positives)
		                           and (chr(n) not in neutrals)])

	stackLeft       = [None]	# stack with each entry containing five
	stackRight      = [None]	# .. elements;  note that entry zero is not
	stackLeftScore  = [None]	# .. used
	stackRightScore = [None]
	stackLower      = [None]
	top   = 0
	score = 0

	for ix in range(len(seq)):
		ch = seq[ix]
		if (ch in blockers):
			# emit intervals

			for sp in range(1,top+1):
				left  = stackLeft [sp] + 1
				right = stackRight[sp]

				while (left < right) and (seq[left]  not in positives): left  += 1
				while (right > left) and (seq[right] not in positives): right -= 1

				right += 1
				if (minLength == None) or (right - left >= minLength):
					yield (left,right)

			#empty stack

			stackLeft       = [None]
			stackRight      = [None]
			stackLeftScore  = [None]
			stackRightScore = [None]
			stackLower      = [None]
			top   = 0
			score = 0
			continue

		if   (ch in positives): weight = 1-rate
		elif (ch in neutrals):  weight = -rate
		else: raise ValueError

		score += weight
		#if ("algorithm" in debug):
		#	print >>sys.stderr, "%3d: %c %5.2f" % (ix, ch, score),

		if (weight < 0):
			#if ("algorithm" in debug):
			#	print >>sys.stderr
			continue

		if (top > 0) and (stackRight[top] == ix-1):
			# add this site to the interval on top of the stack

			stackRight     [top] = ix
			stackRightScore[top] = score

			#if ("algorithm" in debug):
			#	print >>sys.stderr, \
			#	      " extending [%d] %d-%d %4.1f %4.1f" \
			#	    % (top,
			#	       stackLeft     [top], stackRight     [top],
			#	       stackLeftScore[top], stackRightScore[top]),

		else:
			# create a one site interval

			top += 1
			if (top >= len(stackLeft)):
				stackLeft       += [None]
				stackRight      += [None]
				stackLeftScore  += [None]
				stackRightScore += [None]
				stackLower      += [None]

			stackLeft      [top] = ix - 1
			stackLeftScore [top] = score - weight
			stackRight     [top] = ix
			stackRightScore[top] = score
			stackLower     [top] = top - 1

			while (stackLower[top] > 0) \
			  and (stackLeftScore[stackLower[top]] > stackLeftScore[top]):
				stackLower[top] = stackLower[stackLower[top]]

			#if ("algorithm" in debug):
			#	print >>sys.stderr, \
			#	      " creating  [%d] %d-%d %4.1f %4.1f -> %d" \
			#	    % (top,
			#	       stackLeft     [top], stackRight     [top],
			#	       stackLeftScore[top], stackRightScore[top],
			#	       stackLower    [top]),

		# merge intervals;  if there is a previous interval with a no-higher
		# left score and no-higher right score, merge this interval (and all
		# intervening ones) into that one

		while (top > 1) \
		  and (stackLower[top] > 0) \
		  and (stackRightScore[stackLower[top]] <= stackRightScore[top]):
			stackRight     [stackLower[top]] = stackRight     [top]
			stackRightScore[stackLower[top]] = stackRightScore[top]
			top = stackLower[top]

			#if ("algorithm" in debug):
			#	print >>sys.stderr, \
			#	      "\n%*s merging   [%d] %d-%d %4.1f %4.1f" \
			#	    % (13, "", top,
			#	       stackLeft[top],      stackRight     [top],
			#	       stackLeftScore[top], stackRightScore[top]),

		#if ("algorithm" in debug):
		#	print >>sys.stderr

	# emit intervals

	for sp in range(1,top+1):
		left  = stackLeft [sp] + 1
		right = stackRight[sp]

		while (left < right) and (seq[left]  not in positives): left  += 1
		while (right > left) and (seq[right] not in positives): right -= 1

		right += 1
		if (minLength == None) or (right - left >= minLength):
			yield (left,right)


if __name__ == "__main__": main()

