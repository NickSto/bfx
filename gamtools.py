#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import json
import errno
import argparse

OUTFMT = {'n':'name', 'i':'identity', 'm':'mapping_quality', 's':'score', 'S':'sequence',
          'q':'quality'}

ARG_DEFAULTS = {'input':sys.stdin, 'outfmt':'im'}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('input', metavar='align.gam.txt', type=argparse.FileType('r'), nargs='?',
    help='The input alignment. Must be in JSON format. Omit to read from stdin.')
  parser.add_argument('-o', '--outfmt',
    help='The stats to print, as a string of letters: "n": name, "i": identity, "m": '
         'mapping_quality, "s": score, "S": sequence, "q": PHRED qualities. Example: "imS" for 3 '
         'columns: identity, mapping_quality, then sequence.')
  parser.add_argument('-j', '--json', action='store_true',
    help='Print the original JSON records instead of the custom, tabular output format.')
  parser.add_argument('-P', '--pretty', action='store_true',
    help='Print a prettified version of the JSON, with spacing and formatting.')
  parser.add_argument('-n', '--names', nargs='+',
    help='Print stats just for the read with these names.')
  parser.add_argument('--names-file', type=argparse.FileType('r'),
    help='Print stats for just the reads with the names in this file. Format: one line per name.')
  parser.add_argument('-i', '--identity', type=float,
    help='Print only alignments with at least this identity.')
  parser.add_argument('-m', '--mapq', type=int,
    help='Print only alignments with at least this mapping_quality.')
  parser.add_argument('-s', '--score', type=int,
    help='Print only alignments with at least this score.')
  parser.add_argument('-N', '--no-nulls', action='store_true',
    help='Omit alignments where one of the output fields is missing, instead of printing "."')

  args = parser.parse_args(argv[1:])

  out_fields = []
  for char in args.outfmt:
    try:
      key = OUTFMT[char]
    except KeyError:
      fail('Error: "{}" not a recognized --outfmt code.'.format(char))
    out_fields.append(key)

  if args.names_file:
    names = read_names(args.names_file)
  elif args.names:
    names = set(args.names)
  else:
    names = set()

  line_num = 0
  for line in args.input:
    line_num += 1
    record = json.loads(line)
    if names and record['name'] not in names:
      continue
    if args.identity and not ('identity' in record and record['identity'] >= args.identity):
      continue
    if args.mapq and not ('mapping_quality' in record and record['mapping_quality'] >= args.mapping_quality):
      continue
    if args.score and not ('score' in record and record['score'] >= args.score):
      continue
    if args.json:
      if args.pretty:
        print(json.dumps(record, indent=2))
      else:
        sys.stdout.write(line)
      continue
    stats = []
    has_null = False
    for key in out_fields:
      if key not in record:
        has_null = True
      stats.append(record.get(key, '.'))
    if args.no_nulls and has_null:
      continue
    print(*stats, sep='\t')


def read_names(names_file):
  names = set()
  for line in names_file:
    name = line.strip()
    names.add(name)
  return names


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except KeyboardInterrupt:
    pass
  except IOError as ioe:
    if ioe.errno != errno.EPIPE:
      raise
