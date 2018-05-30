#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import argparse
import collections
import getreads

QUANT_ORDER = 5
USAGE = "%(prog)s [options] [input_1.fq [input_2.fq output_1.fq output_2.fq]]"
DESCRIPTION = """Trim the 5' ends of reads by sequence content, e.g. by GC content or presence of
N's."""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION, usage=USAGE)
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
  parser.add_argument('-b', '--filt-bases', default='N',
    help='The bases to filter on. Case-insensitive. Default: %(default)s.')
  parser.add_argument('-t', '--thres', type=float, default=0.5,
    help='The threshold. The read will be trimmed once the proportion of filter bases in the '
         'window exceed this fraction (not a percentage). Default: %(default)s.')
  parser.add_argument('-w', '--window', dest='win_len', type=int, default=1,
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
  parser.add_argument('-q', '--quiet', action='store_true',
    help='Don\'t print trimming stats on completion.')
  parser.add_argument('-T', '--tsv', dest='stats_format', default='human',
                      action='store_const', const='tsv')
  return parser


def main(argv):
  parser = make_argparser()
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
  filters = {'win_len':args.win_len, 'thres':args.thres, 'filt_bases':filt_bases, 'invert':invert,
             'min_len':args.min_length}
  try:
    stats = trim_reads(file1_parser, file2_parser, args.outfile1, args.outfile2,
                       filetype1, filetype2, paired, filters)
  finally:
    for filehandle in (args.infile1, args.infile2, args.outfile1, args.outfile2):
      if filehandle and filehandle is not sys.stdin and filehandle is not sys.stdout:
        filehandle.close()

  if not args.quiet:
    print_stats(stats, args.stats_format)


def trim_reads(file1_parser, file2_parser, outfile1, outfile2, filetype1, filetype2, paired,
               filters):
  """Trim all the reads in the input file(s), writing to the output file(s)."""
  min_len = filters['min_len']
  trims1 = collections.Counter()
  trims2 = collections.Counter()
  omitted1 = collections.Counter()
  omitted2 = collections.Counter()
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
    read1, trim_len1 = trim_read(read1, filters, filetype1)
    trims1[trim_len1] += 1
    if paired:
      read2, trim_len2 = trim_read(read2, filters, filetype2)
      trims2[trim_len2] += 1
      # Output reads if they both pass the minimum length threshold (if any was given).
      if min_len is None or (len(read1.seq) >= min_len and len(read2.seq) >= min_len):
        write_read(outfile1, read1, filetype1)
        write_read(outfile2, read2, filetype2)
      else:
        if len(read1.seq) < min_len:
          omitted1[trim_len1] += 1
        if len(read2.seq) < min_len:
          omitted2[trim_len2] += 1
    else:
      # Output read if it passes the minimum length threshold (if any was given).
      if min_len is None or len(read1.seq) >= min_len:
        write_read(outfile1, read1, filetype1)
      else:
        omitted1[trim_len1] += 1
  # Compile stats.
  stats = {}
  stats['reads'] = sum(trims1.values()) + sum(trims2.values())
  stats['trimmed'] = stats['reads'] - trims1[0] - trims2[0]
  stats['omitted'] = sum(omitted1.values()) + sum(omitted2.values())
  if paired:
    stats['trimmed1'] = stats['reads']//2 - trims1[0]
    stats['trimmed2'] = stats['reads']//2 - trims2[0]
    stats['omitted1'] = sum(omitted1.values())
    stats['omitted2'] = sum(omitted2.values())
  # Quintiles for trim lengths.
  stats['quants'] = {'order':QUANT_ORDER}
  if paired:
    stats['quants']['trim1'] = get_counter_quantiles(trims1, order=QUANT_ORDER)
    stats['quants']['trim2'] = get_counter_quantiles(trims2, order=QUANT_ORDER)
    stats['quants']['trim'] = get_counter_quantiles(trims1 + trims2, order=QUANT_ORDER)
    stats['quants']['omitted_trim1'] = get_counter_quantiles(omitted1, order=QUANT_ORDER)
    stats['quants']['omitted_trim2'] = get_counter_quantiles(omitted2, order=QUANT_ORDER)
    stats['quants']['omitted_trim'] = get_counter_quantiles(omitted1 + omitted2, order=QUANT_ORDER)
  else:
    stats['quants']['trim'] = get_counter_quantiles(trims1)
    stats['quants']['omitted_trim'] = get_counter_quantiles(omitted1)
  return stats


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


def trim_read(read, filters, filetype):
  trimmed_seq = trim_seq(read.seq, **filters)
  trim_len = len(read.seq) - len(trimmed_seq)
  read.seq = trimmed_seq
  if filetype == 'fastq':
    # If the output filetype is FASTQ, trim the quality scores too.
    # If there are no input quality scores (i.e. the input is FASTA), use dummy scores instead.
    # "z" is the highest alphanumeric score (PHRED 89), higher than any expected real score.
    qual = read.qual or 'z' * len(read.seq)
    read.qual = qual[:len(read.seq)]
  return read, trim_len


def trim_seq(seq, win_len=1, thres=1.0, filt_bases='N', invert=False, **kwargs):
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
    # Are we over the threshold?
    if bad_bases_count > max_bad_bases:
      break
  # If we exceeded the threshold, trim the sequence at the first (leftmost) bad base in the window.
  if bad_bases_count > max_bad_bases:
    first_bad_base = bad_bases_coords[0]
    return seq[0:first_bad_base]
  else:
    return seq


def get_counter_quantiles(counter, order=5):
  """Return an arbitrary set of quantiles (including min and max values).
  `counter` is a collections.Counter.
  `order` is which quantile to perform (4 = quartiles, 5 = quintiles).
  Warning: This expects a counter which has counted at least `order` elements.
  If it receives a counter with fewer elements, it will simply return `list(counter.elements())`.
  This will have fewer than the usual order+1 elements, and may not fit normal expectations of
  what "quantiles" should be."""
  quantiles = []
  total = sum(counter.values())
  if total <= order:
    return list(counter.elements())
  span_size = total / order
  # Sort the items and go through them, looking for the one at the break points.
  items = list(sorted(counter.items(), key=lambda i: i[0]))
  quantiles.append(items[0][0])
  total_seen = 0
  current_span = 1
  cut_point = int(round(current_span*span_size))
  for item, count in items:
    total_seen += count
    if total_seen >= cut_point:
      quantiles.append(item)
      current_span += 1
      cut_point = int(round(current_span*span_size))
  return quantiles


def print_stats(stats, format='human'):
  if format == 'human':
    lines = get_stats_lines_human(stats)
  elif format == 'tsv':
    lines = get_stats_lines_tsv(stats)
  else:
    fail('Error: Unrecognized format {!r}'.format(format))
  sys.stderr.write('\n'.join(lines).format(**stats)+'\n')


def get_stats_lines_human(stats):
  # Single-stat lines:
  lines = [
    'Total reads in input:\t{reads}',
    'Reads trimmed:\t{trimmed}'
  ]
  if 'trimmed1' in stats and 'trimmed2' in stats:
    lines.append('  For mate 1:\t{trimmed1}')
    lines.append('  For mate 2:\t{trimmed2}')
  lines.append('Reads filtered out:\t{omitted}')
  if 'omitted1' in stats and 'omitted2' in stats:
    lines.append('  For mate 1:\t{omitted1}')
    lines.append('  For mate 2:\t{omitted2}')
  # Quantile lines:
  quantile_lines = [
    ('Bases trimmed quintiles', 'trim'),
    ('  For mate 1', 'trim1'),
    ('  For mate 2', 'trim2'),
    ('Bases trimmed quantiles from filtered reads', 'omitted_trim'),
    ('  For mate 1', 'omitted_trim1'),
    ('  For mate 2', 'omitted_trim2')
  ]
  for desc, stat_name in quantile_lines:
    if stat_name in stats['quants']:
      quants_values = stats['quants'][stat_name]
      if quants_values:
        quants_str = ', '.join(map(str, quants_values))
      else:
        quants_str = 'N/A'
      line = desc+':\t'+quants_str
      lines.append(line)
  return lines


def get_stats_lines_tsv(stats):
  lines = ['{reads}']
  if 'trimmed1' in stats and 'trimmed2' in stats:
    lines.append('{trimmed}\t{trimmed1}\t{trimmed2}')
  else:
    lines.append('{trimmed}')
  if 'omitted1' in stats and 'omitted2' in stats:
    lines.append('{omitted}\t{omitted1}\t{omitted2}')
  else:
    lines.append('{omitted}')
  for stat_name in ('trim', 'trim1', 'trim2', 'omitted_trim', 'omitted_trim1', 'omitted_trim2'):
    if stat_name in stats['quants']:
      quants_values = stats['quants'][stat_name]
      lines.append('\t'.join(map(str, quants_values)))
  return lines


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
