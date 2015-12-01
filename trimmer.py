#!/usr/bin/env python
from __future__ import division
import sys
import argparse
import getreads

OPT_DEFAULTS = {'win_len':10, 'thres':1.0, 'filt_bases':'N'}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Filter whole reads or trim reads by various criteria like quality scores or N
content."""

def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
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
  parser.add_argument('-A', '--acgt', action='store_true',
    help='Filter on any non-ACGT base (shortcut for "--invert --filt-bases ACGT").')
  parser.add_argument('-I', '--iupac', action='store_true',
    help='Filter on any non-IUPAC base (shortcut for "--invert --filt-bases ACGTUWSMKRYBDHVN-").')

  args = parser.parse_args(argv[1:])

  if args.infile1 and args.infile2 and not (args.outfile1 and args.outfile2):
    fail('Error: If giving two input files (paired end), must specify both output files.')

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

  if args.out_filetype:
    filetype1 = args.out_filetype
    filetype2 = args.out_filetype

  filt_bases = args.filt_bases
  invert = args.invert
  if args.acgt:
    filt_bases = 'ACGT'
    invert = True
  elif args.iupac:
    filt_bases = 'ACGTUWSMKRYBDHVN-'
    invert = True

  trim_reads(file1_parser, file2_parser, args.outfile1, args.outfile2, filetype1, filetype2, paired,
             args.win_len, args.thres, filt_bases, invert)

  if args.infile1 and args.infile1 is not sys.stdin:
    args.infile1.close()
  if paired:
    args.infile2.close()


def trim_reads(file1_parser, file2_parser, outfile1, outfile2, filetype1, filetype2, paired,
               win_len, thres, filt_bases, invert):
  read1 = None
  read2 = None
  while True:
    # Read in reads.
    try:
      read1 = next(file1_parser)
      if paired:
        read2 = next(file2_parser)
    except StopIteration:
      break
    # Do trimming.
    read1.seq = trim(read1.seq, win_len, thres, filt_bases, invert)
    read1.qual = read1.qual[:len(read1.seq)]
    if paired:
      read2.seq = trim(read2.seq, win_len, thres, filt_bases, invert)
      read2.qual = read2.qual[:len(read2.seq)]
      # Output reads if they both passed the filters.
      if read1.seq and read2.seq:
        write_read(outfile1, read1, filetype1)
        write_read(outfile2, read2, filetype2)
    else:
      # Output read if it passed the filters.
      if read1.seq:
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


def trim(seq, win_len, thres, filt_bases, invert):
  max_bad_bases = win_len * thres
  window = []
  bad_bases_count = 0
  bad_bases_coords = []
  for coord, base in enumerate(seq.upper()):
    # Shift window, adjust bad_bases_count and bad_bases_coords list.
    window.append(base)
    # Is the new base a bad base?
    if invert:
      bad_base = base not in filt_bases
    else:
      bad_base = base in filt_bases
    if bad_base:
      bad_bases_count += 1
      bad_bases_coords.append(coord)
    if len(window) > win_len:
      first_base = window.pop(0)
      # Is the first base in the window a bad base?
      if invert:
        bad_base = first_base not in filt_bases
      else:
        bad_base = first_base in filt_bases
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
