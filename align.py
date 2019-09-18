#!/usr/bin/env python3
import argparse
import distutils.spawn
import distutils.version
import logging
import os
import pathlib
import subprocess
import sys
assert sys.version_info.major >= 3, 'Python 3 required'

ALIGNER_DATA = {
  'bowtie2': {
    'required': ('bowtie2-build', 'bowtie2'),
    'index-endings': ('1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2'),
    'opts': [],
  },
  'bwa': {
    'required': ('bwa',),
    'index-endings': ('amb', 'ann', 'bwt', 'pac', 'sa'),
    'opts': ['-M'],
  },
}
DESCRIPTION = """Align reads.
This will automatically run all the commands needed to get you from raw reads to an indexed BAM
file, including indexing the reference sequence."""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('aligner', choices=ALIGNER_DATA.keys(),
    help="Aligner to use. 'bwa' will use BWA-MEM.")
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
  parser.add_argument('-R', '--ref-base',
    help='The base path for the reference index. E.g. "-R path/to/ref" if your index files are '
      'named like "path/to/ref.bwt" and "path/to/ref.pac". Default: the same path as the actual '
      'reference file.')
  parser.add_argument('-i', '--keep-index', action='store_true',
    help="Don't delete the reference index files. Default: leave it as you found it "
      "(clean up if they didn't exist, leave them if they did).")
  parser.add_argument('-I', '--delete-index', action='store_true',
    help='Delete the reference index files, even if they already existed.')
  parser.add_argument('-t', '--threads', type=int, default=1,
    help='Threads to use when aligning. For bowtie2, this also speeds up indexing. '
      'Default: %(default)s')
  opts_str = ''
  for aligner, data in ALIGNER_DATA.items():
    opts = data.get('opts')
    if opts:
      if opts_str:
        opts_str += ', '
      opts_str += "'"+' '.join(opts)+f"' for '{aligner}'"
  parser.add_argument('-O', '--aligner-opts', type=split_opt_list,
    help="Options to pass to the alignment command. If giving a single option, prepend '|' to "
      "keep it from being parsed as an argument to this script. E.g. \"-O '|-M'\". "
      'Defaults: '+opts_str+', none for the rest.')
  parser.add_argument('-N', '--name-sort', dest='sort_key', action='store_const', const='name',
    default='coord',
    help='Sort the output BAM by name instead of coordinate.')
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

  # Check for required commands.
  missing = []
  for command in ALIGNER_DATA[args.aligner]['required']:
    if not distutils.spawn.find_executable(command):
      missing.append(command)
  if missing:
    fail("Error: Missing required command(s) '"+"', '".join(missing)+"'.")

  # Determine paths.
  if args.ref_base:
    ref_base = args.ref_base
  else:
    ref_base = str(args.ref)
  format = get_format(args.out, args.format)
  sam_path, out_path = get_paths(args.reads1, args.out, format)
  if not args.clobber:
    for path in sam_path, out_path:
      if path.exists():
        fail(f'Error: output path {path} already exists. Use -c/--clobber to overwrite.')

  was_indexed = is_indexed(args.aligner, ref_base)

  # Index reference (if needed).
  if not was_indexed:
    clear_indices(args.aligner, ref_base)
    index_ref(args.aligner, args.ref, ref_base, threads=args.threads)

  # Do alignment.
  align(
    args.aligner, ref_base, args.reads1, args.reads2, sam_path, threads=args.threads,
    opts=args.aligner_opts
  )

  if format == 'bam':
    # Convert output.
    # $ samtools view -Sb align.sam | samtools sort -o - dummy > align.bam
    convert(sam_path, out_path, sort_key=args.sort_key)

    # Index output.
    # $ samtools index align.bam
    if args.sort_key == 'coord':
      index_bam(out_path)

    # Clean up intermediate file(s).
    os.remove(sam_path)

  if (not was_indexed and not args.keep_index) or args.delete_index:
    clear_indices(args.aligner, ref_base)

  if out_path.is_file():
    logging.error(f'Success! Output is in {str(out_path)!r}')
  else:
    fail(f'Error: Output file missing {str(out_path)!r}')


def get_paths(reads1_arg, out_arg, format):
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
  return sam_path, out_path


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


def is_indexed(aligner, ref_base):
  for ending in ALIGNER_DATA[aligner]['index-endings']:
    ref_path_str = ref_base+'.'+ending
    if os.path.exists(ref_path_str):
      if not os.path.isfile(ref_path_str):
        fail(f'Error: Index path {ref_path_str!r} exists, but is not a regular file.')
    else:
      return False
  return True


def clear_indices(aligner, ref_base):
  for ending in ALIGNER_DATA[aligner]['index-endings']:
    ref_path_str = ref_base+'.'+ending
    if os.path.isfile(ref_path_str):
      os.remove(ref_path_str)


def index_ref(aligner, ref, ref_base, threads=1):
  kwargs = {}
  if aligner == 'bowtie2':
    # $ bowtie2-build ref.fa ref
    cmd = ['bowtie2-build']
    version = get_bowtie2_version()
    if version is None:
      logging.warning('Warning: Unable to determine bowtie2 version.')
    elif version >= distutils.version.LooseVersion('2.2.7'):
      cmd.extend(['--threads', threads])
    cmd.extend([ref, ref_base])
    kwargs['stdout'] = subprocess.DEVNULL
  elif aligner == 'bwa':
    # $ bwa index -a algo -p ref ref.fa
    # Use the 'is' indexing algorithm, unless the reference is > 2GB.
    if os.path.getsize(ref) < 2000000000:
      algorithm = 'is'
    else:
      algorithm = 'bwtsw'
    cmd = ('bwa', 'index', '-a', algorithm, '-p', ref_base, ref)
  run_command(cmd, **kwargs)


def align(aligner, ref_base, reads1_path, reads2_path, sam_path, threads=1, opts=None):
  if opts is None:
    opts = ALIGNER_DATA[aligner]['opts']
  if aligner == 'bowtie2':
    # $ bowtie2 -x ref -1 reads_1.fq -2 reads_2.fq -S align.sam
    cmd = (
      ['bowtie2'] + opts +
      ['-p', threads, '-x', ref_base, '-1', reads1_path, '-2', reads2_path, '-S', sam_path]
    )
    run_command(cmd)
  elif aligner == 'bwa':
    # $ bwa mem [opts] ref reads_1.fq reads_2.fq > align.sam
    cmd = ['bwa', 'mem'] + opts + ['-t', threads, ref_base, reads1_path, reads2_path]
    with sam_path.open('wb') as sam_file:
      run_command(cmd, stdout=sam_file)


def convert(sam_path, bam_path, sort_key='coord'):
  cmd1 = ['samtools', 'view', '-Sb', sam_path]
  cmd2 = ['samtools', 'sort', '-o', '-', 'dummy']
  if sort_key == 'name':
    cmd2[2:2] = ['-n']
  logging.warning(
    '$ '+' '.join(map(str, cmd1))+' \\\n'
    +'  | '+' '.join(map(str, cmd2))+' \\\n'
    +'  > '+str(bam_path)
  )
  proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
  with bam_path.open('wb') as bam_file:
    proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=bam_file)
  proc1.stdout.close()
  for cmd, proc in zip((cmd1, cmd2), (proc1, proc2)):
    if proc.wait():
      fail(f'Error: {cmd[0]} {cmd[1]} returned with a non-zero exit code.')


def index_bam(bam_path):
  index_path = bam_path.parent / (bam_path.name+'.bai')
  if index_path.exists():
    os.remove(index_path)
  cmd = ('samtools', 'index', bam_path)
  if bam_path.is_file():
    run_command(cmd)
  else:
    fail(f'Error: Output BAM missing ({bam_path})')


def get_samtools_version(exe='samtools'):
  cmd = (exe,)
  result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
  if result.returncode != 1:
    return None
  output = str(result.stderr, 'utf8')
  # Find the version line.
  line_fields = None
  for line in output.splitlines():
    if line.lower().startswith('version:'):
      line_fields = line.split()
  if line_fields is None:
    return None
  # Find the version number in the line.
  ver_str = line_fields[1]
  # Verify it looks like a version number: does it start with two decimal-separated integers?
  ver_fields = ver_str.split('.')
  if len(ver_fields) <= 1:
    return None
  try:
    int(ver_fields[0])
    int(ver_fields[1])
  except ValueError:
    return None
  logging.info(f'Info: Successfully determined samtools version to be {ver_str}.')
  return distutils.version.LooseVersion(ver_str)


def get_bowtie2_version(exe='bowtie2-build'):
  cmd = (exe, '--version')
  result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
  if result.returncode != 0:
    return None
  output = str(result.stdout, 'utf8')
  # Find the version line.
  line_fields = None
  for line_num, line in enumerate(output.splitlines()):
    if line_num == 0:
      line_fields = line.split()
      if not line_fields:
        continue
      exe_path = pathlib.Path(line_fields[0])
      if exe_path.name.startswith('bowtie2-build'):
        break
  if line_fields is None:
    return None
  # Find the version number in the line.
  ver_str = line_fields[-1]
  # Verify it looks like the right version number:
  # Is it at least two decimal-delimited fields, the first of which is 2?
  ver_fields = ver_str.split('.')
  if len(ver_fields) > 1 and ver_fields[0] == '2':
    logging.info(f'Info: Successfully determined bowtie2 version to be {ver_str}.')
    return distutils.version.LooseVersion(ver_str)
  else:
    return None


def run_command(cmd, **kwargs):
  cmd_strs = [str(arg) for arg in cmd]
  logging.warning('$ '+' '.join(cmd_strs))
  if subprocess.call(cmd_strs, **kwargs):
    fail(f'Error: {cmd[0]} returned with a non-zero exit code.')


def split_opt_list(opts_str):
  opts_str2 = opts_str.lstrip('|')
  return opts_str2.split()


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
