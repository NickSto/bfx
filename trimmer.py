#!/usr/bin/env python
from __future__ import division
import sys
import argparse
import getreads

OPT_DEFAULTS = {'win_len':1, 'thres':1.0, 'filt_bases':'N'}
USAGE = "%(prog)s [options] [input_1.fq [input_2.fq output_1.fq output_2.fq]]"
DESCRIPTION = """Trim the 5' ends of reads by sequence content, e.g. by GC content or presence of
N's."""

def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION, usage=USAGE)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile1', metavar='reads_1.fq', nargs='?', type=argparse.FileType('r'),
    default=sys.stdin,
    help='Input reads (mate 1). Omit to read from stdin.')
  parser.add_argument('infile2', metavar='reads_2.fq', nargs='?', type=argparse.FileType('r'),
    help='Input reads (mate 2). If given, it will preserve pairs (if one read is filtered out '
         'entirely, the other will also be lost).')
  parser.add_argument('outfile1', metavar='reads.filt_1.fq', nargs='?', type=argparse.FileType('w'),
    default=sys.stdout,
    help='Output file for mate 1. WARNING: Will overwrite.')
  parser.add_argument('outfile2', metavar='reads.filt_2.fq', nargs='?', type=argparse.FileType('w'),
    help='Output file for mate 2. WARNING: Will overwrite.')
  parser.add_argument('-f', '--format', dest='filetype', choices=('fasta', 'fastq'),
    help='Input read format.')
  parser.add_argument('-F', '--out-format', dest='out_filetype', choices=('fasta', 'fastq'),
    help='Output read format. Default: whatever the input format is.')
  parser.add_argument('-b', '--filt-bases',
    help='The bases to filter on. Case-insensitive. Default: %(default)s.')
  parser.add_argument('-t', '--thres', type=float,
    help='The threshold. The read will be trimmed once the proportion of filter bases in the '
         'window exceed this fraction (not a percentage). Default: %(default)s.')
  parser.add_argument('-w', '--window', dest='win_len', type=int,
    help='Window size for trimming. Default: %(default)s.')
  parser.add_argument('-i', '--invert', action='store_true',
    help='Invert the filter bases: filter on bases NOT present in the --filt-bases.')
  parser.add_argument('-m', '--min-length', type=int,
    help='Set a minimum read length. Reads which are trimmed below this length will be filtered '
         'out (omitted entirely from the output). Read pairs will be preserved: both reads in a '
         'pair must exceed this length to be kept. Set to 0 to only omit empty reads.')
  parser.add_argument('-A', '--acgt', action='store_true',
    help='Filter on any non-ACGT base (shortcut for "--invert --filt-bases ACGT").')
  parser.add_argument('-I', '--iupac', action='store_true',
    help='Filter on any non-IUPAC base (shortcut for "--invert --filt-bases ACGTUWSMKRYBDHVN-").')

  args = parser.parse_args(argv[1:])

  # Catch invalid argument combinations.
  if args.infile1 and args.infile2 and not (args.outfile1 and args.outfile2):
    fail('Error: If giving two input files (paired end), must specify both output files.')
  # Determine filetypes, open input file parsers.
  filetype1 = get_filetype(args.infile1, args.filetype)
  file1_parser = iter(getreads.getparser(args.infile1, filetype=filetype1))
  if args.infile2:
    paired = True
    filetype2 = get_filetype(args.infile2, args.filetype)
    file2_parser = iter(getreads.getparser(args.infile2, filetype=filetype2))
  else:
    filetype2 = None
    file2_parser = None
    paired = False
  # Override output filetypes if it was specified on the command line.
  if args.out_filetype:
    filetype1 = args.out_filetype
    filetype2 = args.out_filetype

  # Determine the filter bases and whether to invert the selection.
  filt_bases = args.filt_bases
  invert = args.invert
  if args.acgt:
    filt_bases = 'ACGT'
    invert = True
  elif args.iupac:
    filt_bases = 'ACGTUWSMKRYBDHVN-'
    invert = True

  # Do the actual trimming.
  try:
    trim_reads(file1_parser, file2_parser, args.outfile1, args.outfile2, filetype1, filetype2,
               paired, args.win_len, args.thres, filt_bases, invert, args.min_length)
  finally:
    for filehandle in (args.infile1, args.infile2, args.outfile1, args.outfile2):
      if filehandle and filehandle is not sys.stdin and filehandle is not sys.stdout:
        filehandle.close()


def trim_reads(file1_parser, file2_parser, outfile1, outfile2, filetype1, filetype2, paired,
               win_len, thres, filt_bases, invert, min_length):
  """Trim all the reads in the input file(s), writing to the output file(s)."""
  read1 = None
  read2 = None
  while True:
    # Read in the reads.
    try:
      read1 = next(file1_parser)
      if paired:
        read2 = next(file2_parser)
    except StopIteration:
      break
    # Do trimming.
    read1.seq = trim_read(read1.seq, win_len, thres, filt_bases, invert)
    if filetype1 == 'fastq':
      # If the output filetype is FASTQ, trim the quality scores too.
      # If there are no input quality scores (i.e. the input is FASTA), use dummy scores instead.
      # "z" is the highest alphanumeric score (PHRED 89), higher than any expected real score.
      qual1 = read1.qual or 'z' * len(read1.seq)
      read1.qual = qual1[:len(read1.seq)]
    if paired:
      read2.seq = trim_read(read2.seq, win_len, thres, filt_bases, invert)
      if filetype2 == 'fastq':
        qual2 = read2.qual or 'z' * len(read2.seq)
        read2.qual = qual2[:len(read2.seq)]
      # Output reads if they both pass the minimum length threshold (if any was given).
      if min_length is None or (len(read1.seq) >= min_length and len(read2.seq) >= min_length):
        write_read(outfile1, read1, filetype1)
        write_read(outfile2, read2, filetype2)
    else:
      # Output read if it passes the minimum length threshold (if any was given).
      if min_length is None or len(read1.seq) >= min_length:
        write_read(outfile1, read1, filetype1)


def get_filetype(infile, filetype_arg):
  if infile is sys.stdin:
    if filetype_arg:
      filetype = filetype_arg
    else:
      fail('Error: You must specify the --format if reading from stdin.')
  elif infile:
    if filetype_arg:
      filetype = filetype_arg
    else:
      if infile.name.endswith('.fa') or infile.name.endswith('.fasta'):
        filetype = 'fasta'
      elif infile.name.endswith('.fq') or infile.name.endswith('.fastq'):
        filetype = 'fastq'
      else:
        fail('Error: Unrecognized file ending on "{}". Please specify the --format.'.format(infile))
  else:
    fail('Error: infile is {}'.format(infile))
  return filetype


def write_read(filehandle, read, filetype):
  if filetype == 'fasta':
    filehandle.write('>{name}\n{seq}\n'.format(**vars(read)))
  elif filetype == 'fastq':
    filehandle.write('@{name}\n{seq}\n+\n{qual}\n'.format(**vars(read)))


def trim_read(seq, win_len, thres, filt_bases, invert):
  """Trim an individual read and return its trimmed sequence.
  This will track the frequency of bad bases in a window of length win_len, and trim once the
  frequency goes below thres. The trim point will be just before the first (leftmost) bad base in
  the window (the first window with a frequency below thres). The "bad" bases are the ones in
  filt_bases if invert is False, or any base NOT in filt_bases if invert is True."""
  # Algorithm:
  # The window is a list which acts as a FIFO. As we scan from the left (3') end to the right (5')
  # end, we append new bases to the right end of the window and pop them from the left end.
  # Each base is only examined twice: when it enters the window and when it leaves it.
  # We keep a running total of the number of bad bases in bad_bases_count, incrementing it when bad
  # bases enter the window and decrementing it when they leave.
  # We also track the location of bad bases in the window with bad_bases_coords so we can figure out
  # where to cut if we have to trim.
  max_bad_bases = win_len * thres
  window = []
  bad_bases_count = 0
  bad_bases_coords = []
  for coord, base in enumerate(seq.upper()):
    # Shift window, adjust bad_bases_count and bad_bases_coords list.
    window.append(base)
    # Is the new base we're adding to the window a bad base?
    if invert:
      bad_base = base not in filt_bases
    else:
      bad_base = base in filt_bases
    # If so, increment the total and add its coordinate to the window.
    if bad_base:
      bad_bases_count += 1
      bad_bases_coords.append(coord)
    if len(window) > win_len:
      first_base = window.pop(0)
      # Is the base we're removing (the first base in the window) a bad base?
      if invert:
        bad_base = first_base not in filt_bases
      else:
        bad_base = first_base in filt_bases
      # If so, decrement the total and remove its coordinate from the window.
      if bad_base:
        bad_bases_count -= 1
        bad_bases_coords.pop(0)
    # print bad_bases_coords
    # Are we over the threshold?
    if bad_bases_count > max_bad_bases:
      break
  # If we exceeded the threshold, trim the sequence at the first (leftmost) bad base in the window.
  if bad_bases_count > max_bad_bases:
    first_bad_base = bad_bases_coords[0]
    return seq[0:first_bad_base]
  else:
    return seq


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
