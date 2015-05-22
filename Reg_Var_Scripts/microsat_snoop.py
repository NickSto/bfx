#!/usr/bin/env python
"""
Snoop thru a fasta file looking for microsatellite repeats of given periods
Output format: length_of_repeat 	 left_flank_length	right_flank_length	repeat_motif	hamming_distance	read_name	read_sequence	read_quality

:Author: Bob Harris (rsharris@bx.psu.edu)
"""
from sys          import argv,stdin,stderr,exit
from string       import maketrans
from md5          import new as md5_new
from pyfracluster import dense_intervals

def usage(s=None):
	message = """
usage: microsat_snoop [fasta_file] [options]
  <fasta_file>                Name of file to read sequences from;  if absent,
                              sequences are read from stdin
  --fast1                     Input file is in fastq format
                              (this is the default)
  --fastq                     Input file is in fastq format
                              (default is fasta unless filename is .fastq)
  --fastq:noquals             Input file is in fastq format, but discard quals
  --period=<length>           (mandatory,cumulative) repeat length(s) to be
                              searched for
                              <length> is expected to be small, less than 10
                              <length> can also be a comma-separated list, or
                              a range <low>..<high>
  --rate=<fraction>           control the candidate repeat interval detector;
                              it will consider intervals with at least
                              <fraction> of matches when shifted by the period;
                              <fraction> is between 0 and 1 and can be either a
                              real number or <n>/<d>
                              (default is 6/7)
  --minlength=<length>        minimum length of intervals reported, in bp
                              (default is 20)
  --progress=<count>          how often to report the sequence we're searching
                              (default is no progress report)
  --allowduplicates           process all input sequences
                              (this is the default)
  --noduplicates              ignore any input sequence that's the same as an
                              earlier sequence
  --nonearduplicates          ignore any input sequence that has the same first
                              100 bp as an earlier sequence
  --nonearduplicate=<length>  ignore any input sequence that has the same first
                              <length> bp as an earlier sequence
  --hamming=<count>           Don't report candidate repeat intervals that have
                              more than <count> mismatches
                              (default is to do no such filtering)
  --prefix=<length>           Don't report candidate repeat intervals that
                              start within <length> of the sequence start
                              (default is to do no such filtering)
  --suffix=<length>           Don't report candidate repeat intervals that
                              end within <length> of the sequence end
                              (default is to do no such filtering)
  --subsample=<k>/<n>         Process only the <k>th sequence of every group of
                              <n> sequences;  <k> ranges from 1 to <n>
  --multipleruns              Consider all candidate intervals in a sequence
                              (default is to consider only the longest)
  --partialmotifs             Consider microatelites with a partial motif
                              (default is to consider only whole motifs)
  --splitbyvalidity           Preprocess sequences, splitting at Ns;  this
                              prevents candidates from including Ns
                              (default is not to split)
  --noflankdisplay            Show entire sequence as flanking regions
                              (this is the default)
  --flankdisplay=<length>     Limit length of flanking regions shown
  --readnamesuffix=<string>   Root of suffix to append to read names;  e.g. 1
                              for forward, 2 for reverse;  this triggers other
                              info to be included in the suffix
                              (default is "1" for fastq;  no suffix for fasta)
  --head=<number>             limit the number of sequences processed
  --markend                   Write a marker line upon completion
                              (default is not to write a marker)
  --help=details              Describe the process, and quit"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


detailedDescription = """In broad terms, the process works as follows:

(1) Identify intervals that are highly correlated with the interval shifted by
    P (the repeat period).  These intervals are called "runs" or "candidates".
    The level of correlation required is controlled by rateThreshold. 
    Depending on whether we want to look for more than one microsat, we either
    find the longest such run (simple algorithm) or many runs (more complicated
    algorithm). The following steps are then performed on each run.

(2) Find the most likely repeat motif in the run.  This is done by counting
    all kmers (of length P) and choosing the most frequent.  If that kmer is
    itself covered by a sub-repeat we discard this run.  The idea is that we
    can ignore a 6-mer like ACGACG because we will find it when we are looking
    for 3-mers.

(3) Once we identify the most likely repeat motif, we then modify the
    interval, adjusting start and end to find the interval that has the fewest
    mismatches vs. a sequence of the motif repeated (hamming distance).  Only
    whole copies of the motif are considered.

(4) At this point we have a valid microsat interval (in the eyes of the
    program). It is subjected to some filtering stages (hamming distance or too
    close to an end), and if it satisfies those conditions, it's reported to
    the user."""

def main():
	global debug

	#=== parse the command line ===

	inputFilename         = None
	inputFormat           = None
	repeatPeriods         = []
	rateThreshold         = 6 / 7.0
	lengthThreshold       = 20
	reportProgress        = None
	discardDuplicates     = False
	discardNearDuplicates = False
	nearDuplicatePrefix   = 100
	hammingThreshold      = None
	prefixThreshold       = None
	suffixThreshold       = None
	subsampleK            = None
	subsampleN            = None
	reportMultipleRuns    = False
	allowPartialMotifs    = False
	splitByValidity       = False
	flankDisplayLimit     = None
	readNameSuffix        = None
	headLimit             = None
	markEndOfFile         = False
	debug                 = []

	for arg in argv[1:]:
		if (arg == "--fasta"):
			inputFormat = "fasta"
		elif (arg == "--fastq"):
			inputFormat = "fastq"
		elif (arg == "--fastq:noquals"):
			inputFormat = "fastq:noquals"
		elif (arg.startswith("--period=")):
			val = arg.split("=",1)[1]
			for period in val.split(","):
				if (".." in period):
					(lowPeriod,highPeriod) = period.split("..",1)
					lowPeriod  = int(lowPeriod)
					highPeriod = int(highPeriod)
					for period in xrange(lowPeriod,highPeriod+1):
						repeatPeriods += [period]
				else:
					repeatPeriods += [int(period)]
		elif (arg.startswith("--rate=")):
			val = arg.split("=",1)[1]
			rateThreshold = float_or_fraction(val)
			assert (0.0 < rateThreshold <= 1.0), "%s not a valid rate" % val
		elif (arg.startswith("--minlength=")):
			val = arg.split("=",1)[1]
			lengthThreshold = int(val)
			assert (lengthThreshold >= 0)
		elif (arg.startswith("--progress=")):
			val = arg.split("=",1)[1]
			reportProgress = int(val)
		elif (arg == "--allowduplicates"):
			discardDuplicates     = False
			discardNearDuplicates = False
		elif (arg == "--noduplicates"):
			discardDuplicates     = True
			discardNearDuplicates = False
		elif (arg == "--nonearduplicates"):
			discardDuplicates     = False
			discardNearDuplicates = True
		elif (arg.startswith("--nonearduplicate=")):
			val = arg.split("=",1)[1]
			discardDuplicates     = False
			discardNearDuplicates = True
			nearDuplicatePrefix   = int(val)
			assert (nearDuplicatePrefix > 0)
		elif (arg.startswith("--hamming=")):
			val = arg.split("=",1)[1]
			hammingThreshold = int(val)
			assert (hammingThreshold >= 0)
		elif (arg.startswith("--prefix=")):
			val = arg.split("=",1)[1]
			prefixThreshold = int(val)
			assert (prefixThreshold >= 0)
		elif (arg.startswith("--suffix=")):
			val = arg.split("=",1)[1]
			suffixThreshold = int(val)
			assert (suffixThreshold >= 0)
		elif (arg.startswith("--subsample=")):
			val = arg.split("=",1)[1]
			(k,n) = val.split("/",2)
			subsampleK = int(k)
			subsampleN = int(n)
			assert (0 < subsampleK <= subsampleN)
		elif (arg == "--multipleruns"):
			reportMultipleRuns = True
		elif (arg == "--partialmotifs"):
			allowPartialMotifs = True
		elif (arg == "--splitbyvalidity"):
			splitByValidity = True
		elif (arg == "--noflankdisplay"):
			flankDisplayLimit = None
		elif (arg.startswith("--flankdisplay=")):
			val = arg.split("=",1)[1]
			flankDisplayLimit = int(val)
			assert (flankDisplayLimit >= 0)
		elif (arg.startswith("--readnamesuffix")):
			readNameSuffix = arg.split("=",1)[1]
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(arg.split("=",1)[1])
		elif (arg == "--markend"):
			markEndOfFile = True
		elif (arg == "--help=details"):
			exit (detailedDescription)
		elif (arg.startswith("--debug=")):
			debug += (arg.split("=",1)[1]).split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (inputFilename == None):
			inputFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	#=== determine periods of interest ===

	if (repeatPeriods == []):
		usage("you gotta give me a repeat period")

	periodSeed = {}
	for period in repeatPeriods:
		if (period < 1): usage("period %d is not valid" % period)
		periodSeed[period] = True

	repeatPeriods = [period for period in periodSeed]
	repeatPeriods.sort()

	#=== determine input format ===

	if   (inputFormat == "fasta"):           sequence_reader = fasta_sequences
	elif (inputFormat == "fastq"):           sequence_reader = fastq_sequences
	elif (inputFormat == "fastq:noquals"):   sequence_reader = fastq_sequences
	elif (inputFilename == None):            sequence_reader = fasta_sequences
	elif (inputFilename.endswith(".fastq")): sequence_reader = fastq_sequences
	elif (inputFilename.endswith(".fq")):    sequence_reader = fastq_sequences
	else:                                    sequence_reader = fasta_sequences

	if (inputFilename != None): inputF = file(inputFilename,"rt")
	else:                       inputF = stdin

	if   (readNameSuffix == None) \
	 and (sequence_reader == fastq_sequences) \
	 and (inputFormat != "fastq:noquals"):
		readNameSuffix = "1"

	#=== process the sequences ===

	sequenceSeen = {}

	numSequences = 0
	for seqInfo in sequence_reader(inputF):
		numSequences += 1
		if (headLimit != None) and (numSequences > headLimit):
			print >>stderr, "limit of %d sequences reached" % headLimit
			break

		if (sequence_reader == fastq_sequences):
			(name,sequence,quals) = seqInfo
			if (inputFormat == "fastq:noquals"): quals = None
		else:
			(name,sequence) = seqInfo
			quals = None

		if (reportProgress != None) and (numSequences % reportProgress == 0):
			print >>stderr, "%s %d" % (name,numSequences)

		# if we're subsampling and not interested in this sequence, skip it

		if (subsampleN != None):
			if ((numSequences-1) % subsampleN != (subsampleK-1)):
				continue

		# if this sequence is shorter than the length of interest, skip it

		seqLen = len(sequence)
		if (seqLen < period) or (seqLen < lengthThreshold): continue

		# if we're not interested in duplicates and this is one, skip it;
		# note that we assume no hash collisions occur, i.e. that all hash
		# matches are truly sequence matches

		if (discardDuplicates):
			h = hash108(sequence)
			if (h in sequenceSeen): continue
			sequenceSeen[h] = True
		elif (discardNearDuplicates):
			h = hash108(sequence[:nearDuplicatePrefix])
			if (h in sequenceSeen): continue
			sequenceSeen[h] = True

		# split the sequence into chunks of valid nucleotides

		if (splitByValidity):
			chunks = [(start,end) for (start,end) in nucleotide_runs(sequence)]
		else:
			chunks = [(0,len(sequence))]

		# evaluate for each period of interest

		for period in repeatPeriods:

			# operate on each chunk

			for (chunkStart,chunkEnd) in chunks:
				chunkLen = chunkEnd - chunkStart
				if (chunkLen < period) or (chunkLen < lengthThreshold): continue

				if ("validity" in debug) or ("correlation" in debug) or ("runs" in debug):
					print >>stderr, ">%s_%d_%d" % (name,chunkStart,chunkEnd)

				# compute correlation sequence

				corr = correlation_sequence(sequence,period,chunkStart,chunkEnd)

				if ("correlation" in debug) or ("runs" in debug):
					print >>stderr, sequence[chunkStart:chunkEnd]
					print >>stderr, corr

				# find runs (candidates for being a microsat) 

				if (reportMultipleRuns):
					runs = all_suitable_runs(corr,lengthThreshold-period,rateThreshold)
				else:
					runs = longest_suitable_run(corr,lengthThreshold,rateThreshold)
				if (runs == []): continue

				if ("runs" in debug):
					for (start,end) in runs:
						run = [" "] * seqLen
						for ix in xrange(start-period,end):
							run[ix] = "*"
						print >>stderr, "".join(run)

				if ("candidates" in debug):
					for (start,end) in runs:
						print >>stderr, "%s %d %d" % (name,start,end)

				# process runs and report those that pass muster

				runCount = 0
				for (start,end) in runs:
					runCount += 1

					start = chunkStart + start - period
					end   = chunkStart + end

					(kmer,d,start,end) = find_repeat_element(period,sequence,start,end,allowPartials=allowPartialMotifs)
					if (kmer == None): continue    # (no useful repeat kmer was found)

					rptExtent = end - start
					prefixLen = start
					suffixLen = seqLen - end
					if (rptExtent <= period): continue
					if (hammingThreshold != None) and (d         > hammingThreshold): continue
					if (prefixThreshold  != None) and (prefixLen < prefixThreshold):  continue
					if (suffixThreshold  != None) and (suffixLen < suffixThreshold):  continue

					if (flankDisplayLimit == None):
						seq = sequence[:start] \
							+ sequence[start:end].lower() \
							+ sequence[end:]
					else:
						seq = sequence[max(chunkStart,start-flankDisplayLimit):start] \
							+ sequence[start:end].lower() \
							+ sequence[end:min(chunkEnd,end+flankDisplayLimit)]
					reportName = name
					if (readNameSuffix != None):
						reportName += "_"+readNameSuffix+"_per"+str(period)+"_"+str(runCount)
					if (quals == None): quals = ""
					else:               quals = "\t" + quals
					print "%d\t%d\t%d\t%s\t%d\t%s\t%s%s" \
						% (rptExtent,prefixLen,suffixLen,kmer,d,reportName,seq,quals)

	if (markEndOfFile):
		print "# microsat_snoop end-of-file"

	if (inputF != stdin):
		inputF.close()


# correlation_sequence--
#	Compute the correlation sequence for a given period.  This is a sequence
#	of + and - indicating whether the base at a given position matches the one
#	P positions earlier (where P is the period).  The first P positions are
#	blank.  Positions with single character runs longer than the period are
#	considered as non-matches, unless the period is 1.

def correlation_sequence(sequence,period,start=None,end=None):
	if (start == None): start = 0
	if (end   == None): end   = len(sequence)

	prevCh = sequence[start]
	run    = 1
	for ix in xrange(start+1,start+period):
		ch = sequence[ix]
		if (ch != prevCh): run =  1
		else:              run += 1
		prevCh = ch

	corr = [" "] * period
	for ix in xrange(start+period,end):
		rptCh = sequence[ix-period]
		ch    = sequence[ix]
		if (ch != prevCh): run =  1
		else:              run += 1
		if    (ch    in "ACGT") \
		  and (ch == rptCh) \
		  and ((period == 1) or (run < period)):
			corr += ["+"]
		else:
			corr += ["-"]
		prevCh = ch

	return "".join(corr)


# longest_suitable_run--
#	Find longest run with a good enough rate (of positives).
#
#	We score a "+" as 1-r and anything else as -r.  This is based on the fol-
#	lowing derivation (p is the number of "+"s, n is the number of non-"+"s):
#		p/(p+n) >= r
#		==> p >= rp + rn
#		==> (1-r)p - rn >= 0
#
#	We adapt an algorithm from "Programming Pearls", pg. 81 (2000 printing).
#
# $$$ we use the denominator as the threshold, but we really should use the
# $$$ .. numerator, comparing it to minLength*rate
#
# $$$ this needs to account for $$$ this situation:
# $$$   sequence: ACGACGACGACGTTATTATTATTA
# $$$   matches:     +++++++++---+++++++++
# $$$ this is currently considered to be one interval (if rate <= 6/7), but it
# $$$ ought to be two;  we can't just post-process, though, because some other
# $$$ interval might be longer than the longest half of this;  maybe what we
# $$$ need to do is consider matches at distances -P and -2P, or if we match
# $$$ -P but that itself was a mismatch, we should carry the mismatch forward

def longest_suitable_run(seq,minLength,rate):

	maxEndingHere = 0
	maxSoFar      = 0
	start         = None

	for ix in xrange(len(seq)):
		if (seq[ix] == "+"): s = 1-rate
		else:                s = -rate

		if (maxEndingHere+s < 0):
			maxEndingHere = 0
			block         = ix
		else:
			maxEndingHere += s
			if (maxEndingHere >= maxSoFar):
				maxSoFar = maxEndingHere
				start    = block + 1
				end      = ix + 1

	if (start == None) or (end - start < minLength):
		return []
	else:
		return [(start,end)]


# all_suitable_runs--
#	Find all non-overlapping runs with a good enough rate (of positives), and
#	which meet our length threshold.
# $$$ this needs to post-process the intervals, splitting them to account for
# $$$ this situation:
# $$$   sequence: ACGACGACGACGTTATTATTATTA
# $$$   matches:     +++++++++---+++++++++
# $$$ this is currently reported as one interval (if rate <= 6/7), but it
# $$$ ought to be two

def all_suitable_runs(seq,minCorrLength,rate):
	return [(start,end) for (start,end) in dense_intervals(seq,rate,"+","-",blockers=None,minLength=minCorrLength)]


# find_repeat_element--
#	Find the most plausible repeat element for a run, and nudge the ends of
#	the run if needed.  Note that we will not consider kmers that represent
#	shorter repeats.  For example, we won't report ACTACT as a 6-mer since we
#	consider this to have a shorter period than 6.

def find_repeat_element(period,seq,start,end,allowPartials=False):

	# count the number of occurences of each k-mer;  note that we can't
	# reject kmers containing smaller repeats yet, since for a sequence like
	# ACACACACACAAACACACACACACACACAC we must first discover ACACAC as the best
	# 6-mer, and THEN reject it;  if we reject ACACAC while counting, we'd end
	# up reporting something like ACACAA as the best motif 

	if ("element" in debug):
		print >>stderr, "find_repeat_element(%d,%d,%d)" % (period,start,end)

	kmerToCount = {}
	kmerToFirst = {}
	for ix in xrange(start,end-(period-1)):
		kmer = seq[ix:ix+period]
		if ("N" in kmer): continue
		if (kmer not in kmerToCount):
			kmerToCount[kmer] = 1
			kmerToFirst[kmer] = ix
		else:
			kmerToCount[kmer] += 1
		#if ("element" in debug):
		#	print >>stderr, "    %d: %s" % (ix,kmer)

	# choose the best k-mer;  this is simply the most frequently occurring one,
	# with ties broken by whichever one came first

	kmers = [(-kmerToCount[kmer],kmerToFirst[kmer],kmer) for kmer in kmerToCount]
	if (kmers == []): return (None,None,start,end)
	kmers.sort()

	if ("element" in debug):
		for (count,first,kmer) in kmers:
			print >>stderr, "    %s: %d" % (kmer,-count)

	(count,first,kmer) = kmers[0]
	if (contains_repeat(kmer)): return (None,None,start,end)

	# determine the hamming distance between the run and a simple repeat, for
	# each "plausible" start and end;  we compute the distance for each such
	# interval, and choose the one with the lowest hamming distance;  ties are
	# broken in a deterministic-but-unspecified manner

	bestD = bestStart = bestEnd = None

	for (s,e) in plausible_intervals(start,end,period,len(seq),allowPartials=allowPartials):
		d = hamming_distance(seq,s,e,kmer)
		if (d == None): continue
		if (bestD == None) or (d <= bestD):
			(bestD,bestStart,bestEnd) = (d,s,e)

	return (kmer,bestD,bestStart,bestEnd)


# plausible_intervals--
#	Yield all plausible intervals intersecting with a run.  We generate all
#	starts within P bp of the run's start.  For each of these, we either (a) try
#	all ends within P bp of run's end, or (b) trim the new interval to a whole
#	multiple of the period, and report this short interval and the longer
#	interval with one more period appended.  Case (a) allows partial motifs,
#	while case (b) only allows whole motifs.

def plausible_intervals(start,end,period,seqLen,allowPartials=False):

	# generate intervals that allow a partial copy of the motif

	if (allowPartials):
		for candStart in xrange(start-(period-1),start+period):
			if (candStart < 0): continue
			for candEnd in xrange(end-(period-1),end+period):
				if (candEnd > seqLen): continue
				if (candEnd <= candStart+period): continue
				yield (candStart,candEnd)

	# -OR- generate intervals that allow only whole copies of the motif

	else:
		for candStart in xrange(start-(period-1),start+period):
			if (candStart < 0): continue
			candEnd = candStart + ((end-candStart)/period)*period
			yield (candStart,candEnd)
			candEnd += period
			if (candEnd <= seqLen): yield (candStart,candEnd)


# hamming_distance--
#	Determine the hamming distance between the run and a simple repeat.
# $$$ improve this by allowing gaps, and stopping when we reach a threshold

kmerToDiffs = {}  # (this is used for memo-ization)

def hamming_distance(seq,start,end,kmer):
	period = len(kmer)
	if (end < start + period): return None

	wholeEnd = start + ((end-start)/period)*period

	if (kmer not in kmerToDiffs):
		kmerToDiffs[kmer] = { kmer:0 }

	d = 0
	for ix in xrange(start,wholeEnd,period):
		qmer = seq[ix:ix+period]	# same size as the kmer motif
		if (qmer in kmerToDiffs[kmer]):
			d += kmerToDiffs[kmer][qmer]
			continue
		diffs = 0
		for iy in xrange(0,period):
			if (qmer[iy] != kmer[iy]): diffs += 1
		kmerToDiffs[kmer][qmer] = diffs
		d += diffs

	if (end > wholeEnd):
		qmer = seq[wholeEnd:end]	# shorter than the kmer motif
		if (qmer in kmerToDiffs[kmer]):
			d += kmerToDiffs[kmer][qmer]
		else:
			diffs = 0
			for iy in xrange(0,len(qmer)):
				if (qmer[iy] != kmer[iy]): diffs += 1
			kmerToDiffs[kmer][qmer] = diffs
			d += diffs

	return d


# fasta_sequences--
#	Read the fasta sequences from a file.  Note that we convert to upper case,
#	and convert any letter other than ACGT to N.

nonDnaMap = maketrans("BDEFHIJKLMOPQRSUVWXYZ","NNNNNNNNNNNNNNNNNNNNN")

def fasta_sequences(f):
	seqName = None
	seqNucs = None

	for line in f:
		line = line.strip()
		if (line.startswith(">")):
			if (seqName != None):
				yield (seqName,"".join(seqNucs))
			seqName = sequence_name(line)
			seqNucs = []
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seqNucs += [line]

	if (seqName != None):
		yield (seqName,"".join(seqNucs).upper().translate(nonDnaMap))


# fastq_sequences--
#	Read the fastq sequences from a file.  Note that we convert to upper case,
#	and convert any letter other than ACGT to N.

def fastq_sequences(f):
	lineNum = 0
	for line in f:
		lineNum += 1
		line = line.strip()

		if (lineNum % 4 == 1):
			assert (line.startswith("@")), \
				   "bad read name at line %d" % lineNum
			seqName = line[1:]
			continue

		if (lineNum % 4 == 2):
			seqNucs = line
			continue

		if (lineNum % 4 == 3):
			assert (line.startswith("+")), \
				   "can't understand line %d:\n%s" % (lineNum,line)
			continue

		quals = line
		assert (len(quals) == len(seqNucs)), \
			   "length mismatch read vs. qualities at line %d" % lineNum
		yield (seqName,"".join(seqNucs).upper().translate(nonDnaMap),quals)

	assert (lineNum % 4 == 0), \
		   "incomplete read at end of file"


# sequence_name--
#	Extract the sequence name from a fasta header.
#	$$$ this may need to be improved $$$

def sequence_name(s):
	s = s[1:].strip()
	if (s == ""): return ""
	else:         return s.split()[0]


# nucleotide_runs--
#	Yield (start,end) for all runs of valid nucleotides in a sequence.

def nucleotide_runs(s):
	runs  = []
	start = None
	for (ix,nuc) in enumerate(s):
		if (nuc in "ACGT"):
			if (start == None):
				start = ix
		else:
			if (start != None):
				yield (start,ix)
				start = None

	if (start != None): yield (start,len(s))


# contains_repeat--
#	Determine whether a short sequence contains a repeated element, such as a
#	6-mer containing a repeated 2-mer (ACACAC) or 3-mer (ACTACT).  The repeat
#	must cover the entire sequence, without mismatches.

def contains_repeat(kmer):
	kmerLength = len(kmer)
	hasRepeat = False
	rptLen = 1
	while (not hasRepeat) and (2 * rptLen <= kmerLength):
		if (kmerLength % rptLen != 0):
			rptLen += 1
			continue
		isRepeat = True
		for i in xrange(rptLen,kmerLength,rptLen):
			if (kmer[i:i+rptLen] != kmer[:rptLen]):
				isRepeat = False
				break
		if (isRepeat):
			hasRepeat = True
			break
		rptLen += 1
	return hasRepeat


# hash108--
#	Return a 108-bit hash "value" of a string

def hash108(s):
	m = md5_new()
	m.update(s)
	return m.hexdigest()[:27]


# float_or_fraction--
#	Convert a string to a number, allowing fractions

def float_or_fraction(s):
	if ("/" in s):
		(numer,denom) = s.split("/",1)
		return float(numer)/float(denom)
	else:
		return float(s)


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return               int(s)   * multiplier
	except ValueError: return int(math.ceil(float(s) * multiplier))


if __name__ == "__main__": main()

