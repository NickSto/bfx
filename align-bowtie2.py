#!/usr/bin/env python3
import argparse
import logging
import os
import pathlib
import subprocess
import sys
assert sys.version_info.major >= 3, 'Python 3 required'

INDEX_ENDINGS = ('1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2')
DESCRIPTION = """Align reads with bowtie2.
This will automatically run all the commands needed to get you from raw reads to an indexed BAM
file, including indexing the reference sequence."""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('ref', metavar='ref.fa', type=pathlib.Path,
    help='Reference sequence.')
  parser.add_argument('reads1', metavar='reads_1.fq', type=pathlib.Path,
    help='Reads, mate 1.')
  parser.add_argument('reads2', metavar='reads_2.fq', type=pathlib.Path,
    help='Reads, mate 2.')
  parser.add_argument('-o', '--out', type=pathlib.Path,
    help='Output file. The aligned BAM will be written to this path.')
  parser.add_argument('-f', '--format', choices=('sam', 'bam'),
    help='Output format.')
  parser.add_argument('-c', '--clobber', action='store_true',
    help='Overwrite intermediate and output files without prompting. Otherwise, the script will '
      'fail if one of these files already exists.')
  parser.add_argument('-i', '--keep-index', action='store_true',
    help="Don't delete the reference index files. Default: leave it as you found it "
      "(clean up if they didn't exist, leave them if they did).")
  parser.add_argument('-I', '--delete-index', action='store_true',
    help='Delete the reference index files, even if they already existed.')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = parser.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  # Determine paths.
  format = get_format(args.out, args.format)
  ref_base, sam_path, out_path = get_paths(args.ref, args.reads1, args.out, format)
  if not args.clobber:
    for path in sam_path, out_path:
      if path.exists():
        fail(f'Error: output path {path} already exists.')

  was_indexed = is_indexed(ref_base)

  # Index reference (if needed).
  # $ bowtie2-build ref.fa ref
  if not was_indexed:
    clear_indices(ref_base)
    index_ref(args.ref, ref_base)

  # Do alignment.
  # $ bowtie2 -x ref -1 reads_1.fq -2 reads_2.fq -S align.sam
  align(ref_base, args.reads1, args.reads2, sam_path)

  if format == 'bam':
    # Convert output.
    # $ samtools view -Sb align.sam | samtools sort -o - dummy > align.bam
    convert(sam_path, out_path)

    # Index output.
    # $ samtools index align.bam
    index_bam(out_path)

    # Clean up intermediate file(s).
    os.remove(sam_path)

  if (not was_indexed and not args.keep_index) or args.delete_index:
    clear_indices(ref_base)

  if out_path.is_file():
    logging.error(f'Success! Output is in {str(out_path)!r}')
  else:
    fail(f'Error: Output file missing {str(out_path)!r}')


def get_paths(ref_arg, reads1_arg, out_arg, format):
  ref_base = get_ref_base(ref_arg)
  base = get_reads_base(reads1_arg)
  if format == 'sam':
    if out_arg:
      sam_path = out_arg
    else:
      sam_path = pathlib.Path(base+'.sam')
    out_path = sam_path
  elif format == 'bam':
    sam_path = pathlib.Path(base+'.sam')
    if out_arg:
      out_path = out_arg
    else:
      out_path = pathlib.Path(base+'.bam')
  return ref_base, sam_path, out_path


def get_ref_base(ref_path):
  basename = ref_path.stem
  return str(ref_path.parent / basename)


def get_reads_base(reads_path):
  basename = reads_path.stem
  if basename.endswith('_1') or basename.endswith('_2'):
    basename = basename[:-2]
  return str(reads_path.parent / basename)


def get_format(out_arg, format_arg):
  if format_arg:
    return format_arg
  if out_arg is not None:
    ext = out_arg.suffix.lower()
    if ext == '.bam':
      return 'bam'
    elif ext == '.sam':
      return 'sam'
  # Default: BAM
  return 'bam'


def is_indexed(ref_base):
  for ending in INDEX_ENDINGS:
    ref_path_str = ref_base+'.'+ending
    if os.path.exists(ref_path_str):
      if not os.path.isfile(ref_path_str):
        fail(f'Error: Index path {ref_path_str!r} exists, but is not a regular file.')
    else:
      return False
  return True


def clear_indices(ref_base):
  for ending in INDEX_ENDINGS:
    ref_path_str = ref_base+'.'+ending
    if os.path.isfile(ref_path_str):
      os.remove(ref_path_str)


def index_ref(ref, ref_base):
  cmd = ('bowtie2-build', ref, ref_base)
  logging.warning('$ '+' '.join(map(str, cmd)))
  subprocess.call(cmd, stdout=subprocess.DEVNULL)


def align(ref_base, reads1_path, reads2_path, sam_path):
  cmd = ('bowtie2', '-x', ref_base, '-1', reads1_path, '-2', reads2_path, '-S', sam_path)
  logging.warning('$ '+' '.join(map(str, cmd)))
  subprocess.call(cmd)


def convert(sam_path, bam_path):
  cmd1 = ('samtools', 'view', '-Sb', sam_path)
  cmd2 = ('samtools', 'sort', '-o', '-', 'dummy')
  logging.warning(
    '$ '+' '.join(map(str, cmd1))+' \\\n'
    +'  | '+' '.join(map(str, cmd2))+' \\\n'
    +'  > '+str(bam_path)
  )
  proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
  with bam_path.open('wb') as bam_file:
    proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=bam_file)
  proc1.stdout.close()


def index_bam(bam_path):
  cmd = ('samtools', 'index', bam_path)
  logging.warning('$ '+' '.join(map(str, cmd)))
  subprocess.call(cmd)


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
