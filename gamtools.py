#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import json
import errno
import argparse
import collections

OUTFMT = collections.OrderedDict((('n', 'name'), ('i', 'identity'), ('m', 'mapping_quality'),
                                  ('s', 'score'), ('S', 'sequence'), ('q', 'quality'),
                                  ('M', 'mappings'), ('e', 'edits'), ('a', 'alts')))
OUTFMT_DESC = OUTFMT.copy()
OUTFMT_DESC['q'] = 'PHRED qualities'
OUTFMT_DESC['M'] = 'number of mappings'
OUTFMT_DESC['e'] = 'number of edits in the top-ranked mapping'
OUTFMT_DESC['a'] = 'alternative sequences of indels'

REVCOMP_TABLE = dict([(ord(a), b) for a, b in zip('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')])

ARG_DEFAULTS = {'input':sys.stdin, 'outfmt':'im'}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Parse, filter, and/or print statistics on alignments in vg's GAM format."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('input', metavar='align.gam.txt', type=argparse.FileType('r'), nargs='?',
    help='The input alignment. Must be in JSON format. Omit to read from stdin.')
  parser.add_argument('-o', '--outfmt',
    help='The stats to print, as a string of letters: '+', '.join(['"{}": {}'.format(char, desc)
          for char, desc in OUTFMT_DESC.items()])+'. Example: "imS" for 3 columns: identity, '
         'mapping_quality, then sequence.')
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
  parser.add_argument('-M', '--mappings', type=int,
    help='Print only alignments with at MOST this many mappings.')
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

  hits = 0
  line_num = 0
  for line in args.input:
    line_num += 1
    record = json.loads(line)
    if names:
      if hits >= len(names):
        break
      elif record['name'] in names:
        hits += 1
      else:
        continue
    if args.identity and not ('identity' in record and record['identity'] >= args.identity):
      continue
    if args.mapq and not ('mapping_quality' in record and record['mapping_quality'] >= args.mapping_quality):
      continue
    if args.score and not ('score' in record and record['score'] >= args.score):
      continue
    if args.mappings and 'path' in record and 'mapping' in record['path'] and len(record['path']['mapping']) > args.mappings:
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
      if key == 'mappings':
        if 'mapping' in record.get('path', {}):
          stats.append(len(record['path']['mapping']))
        else:
          has_null = True
          stats.append('.')
      elif key == 'edits':
        if record.get('path', {}).get('mapping'):
          mapping = get_best_mapping(record['path']['mapping'])
          # We're assuming there's always an 'edit' key in a mapping, since it's too complicated
          # to do otherwise.
          edits = len(mapping['edit'])
          stats.append(edits)
        else:
          has_null = True
          stats.append('.')
      elif key == 'alts':
        if record.get('path', {}).get('mapping'):
          mapping = get_best_mapping(record['path']['mapping'])
          alts = []
          # Same assumption here about the 'edit' key.
          for edit in mapping['edit']:
            # Add any non-SNP variants.
            if edit.get('to_length') != edit.get('from_length') and 'sequence' in edit:
              if mapping.get('position', {}).get('is_reverse'):
                alts.append(get_revcomp(edit['sequence']))
              else:
                alts.append(edit['sequence'])
          stats.extend(alts)
        else:
          has_null = True
          stats.append('.')
      else:
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


def get_best_mapping(mappings):
  best_mapping = None
  best_rank = 0
  for mapping in mappings:
    if mapping['rank'] > best_rank:
      best_mapping = mapping
      best_rank = mapping['rank']
  return best_mapping


def get_revcomp(seq):
  return seq.translate(REVCOMP_TABLE)[::-1]


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
