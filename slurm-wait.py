#!/usr/bin/env python3
import argparse
import getpass
import logging
import pathlib
import subprocess
import sys
import time
assert sys.version_info.major >= 3, 'Python 3 required'

THRESHOLDS = {
  'min_idle_nodes':int,
  'min_idle_cpus':int,
  'min_node_size':int,
  'min_node_size_cpus':int,
  'min_node_size_nodes':int,
  'max_jobs':int,
  'min_jobs':int,
}
USER = getpass.getuser()
PAUSE_TIME = 10
USAGE = """$ %(prog)s [thresholds]
       $ %(prog)s @path/to/args.txt"""
DESCRIPTION = """Determine whether we should keep launching slurm jobs, based on available
resources. Launch this and it will sleep until enough resources are available."""
EPILOG = """You can also give arguments via a file. Just give the path to the file as an argument,
prepended by '@'. The file will be read, and each line will be interpreted as a regular argument to
this script. Remember to put literally every argument on its own line, so '--min-jobs 8' should be
given as two lines: '--min-jobs' and '8'."""


def make_argparser():
  parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG,
    fromfile_prefix_chars='@')
  parser.add_argument('-o', '--output', choices=('min', 'max'),
    help="Choose a node with free resources to run on and print it to stdout. Give either 'min' or "
      "'max' to indicate whether to choose the node with the most or least number of idle CPUs. "
      'Give --cpus-req to indicate how many CPUs the job requires. Only nodes with at least this '
      'many CPUs will be considered.')
  parser.add_argument('-C', '--cpus', type=int, default=1,
    help="How many CPUs are required by the job we're waiting to start. Used when choosing a node "
      'to print with --output.')
  parser.add_argument('-q', '--wait-for-job',
    help="Wait until the job with this name has begun. Useful if you just launched one and don't "
      "want to keep queueing jobs if they're not running.")
  parser.add_argument('-n', '--min-idle-nodes', default=0,
    help='Keep this many nodes idle: if only this many are idle, wait. Can give a number or a file '
      'containing the number (and nothing else).')
  parser.add_argument('-c', '--min-idle-cpus', default=0,
    help='Keep this many CPUs idle: if only this many are idle, wait. Can give a number or a file '
      'containing the number (and nothing else).')
  parser.add_argument('-s', '--min-node-size', default=1,
    help="Minimum node size when counting available resources for the above thresholds. Don't "
      'consider nodes with fewer than this many CPUs. Can give a number or a file containing the '
      'number (and nothing else).')
  parser.add_argument('--min-node-size-cpus', default=0,
    help='Same as --min-node-size, but only for when counting idle cpus for --min-idle-cpus')
  parser.add_argument('--min-node-size-nodes', default=0,
    help='Same as --min-node-size, but only for when counting idle nodes for --min-idle-nodes')
  parser.add_argument('-J', '--max-jobs',
    help="Don't let yourself have more than this many jobs running at once. If you have this many "
      'jobs running, wait. Can give a number or a file containing the number (and nothing else).')
  parser.add_argument('-j', '--min-jobs', default=0,
    help='Always let yourself run at least this many jobs. Even if too few resources are '
      "available and you'd normally wait, keep going if fewer than this many jobs are running. "
      'Note: Does not override --pause-file. Can give a number or a file containing the '
      'number (and nothing else).')
  parser.add_argument('-p', '--pause-file', type=pathlib.Path,
    help='Wait if this file exists. Can be used as a manual pause button.')
  parser.add_argument('-i', '--check-interval', type=int, default=15,
    help='How many seconds to wait between checks for available resources.')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = parser.add_mutually_exclusive_group()
  volume.add_argument('-Q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  thresholds = Thresholds(args)
  thresholds.subdivide_thres('min_node_size', ('min_node_size_cpus', 'min_node_size_nodes'))

  if thresholds.max_jobs is not None and thresholds.min_jobs >= thresholds.max_jobs:
    fail(f'Error: --min-jobs must be < --max-jobs ({thresholds.min_jobs} >= {thresholds.max_jobs}).')

  wait_for_job(args.wait_for_job, args.check_interval)
  wait_for_jobs(thresholds, args.check_interval)
  if count_running_jobs() >= thresholds.min_jobs:
    wait_for_nodes(thresholds, args.check_interval)
    wait_for_cpus(thresholds, args.check_interval)
  wait_for_pause_file(args.pause_file, args.check_interval)

  if args.output:
    node = get_idle('node', args.cpus, args.output)
    print(node)


def wait_for_job(job_name, check_interval):
  if job_name:
    first_loop = True
    printed = False
    while count_running_jobs(job_name) <= 0:
      if not printed and not first_loop:
        print('Waiting for job to begin..')
        printed = True
      if first_loop:
        time.sleep(5)
      else:
        time.sleep(check_interval)
      first_loop = False


def wait_for_jobs(thresholds, check_interval):
  if thresholds.max_jobs is not None:
    printed = False
    while count_running_jobs() >= thresholds.max_jobs:
      if not printed:
        print('Too many jobs running ({} >= {})'.format(count_running_jobs(), thresholds.max_jobs))
        printed = True
      time.sleep(check_interval)


def wait_for_nodes(thresholds, check_interval):
  if thresholds.min_idle_nodes:
    printed = False
    while get_idle('nodes', thresholds.min_node_size_nodes) < thresholds.min_idle_nodes:
      if not printed:
        print(
          'Too few nodes idle ({} < {})'
          .format(get_idle('nodes', thresholds.min_node_size_nodes), thresholds.min_idle_nodes)
        )
        printed = True
      time.sleep(check_interval)


def wait_for_cpus(thresholds, check_interval):
  if thresholds.min_idle_cpus:
    printed = False
    while get_idle('cpus', thresholds.min_node_size_cpus) < thresholds.min_idle_cpus:
      if not printed:
        print(
          'Too few cpus idle ({} < {})'
          .format(get_idle('cpus', thresholds.min_node_size_cpus), thresholds.min_idle_cpus)
        )
        printed = True
      time.sleep(check_interval)


def wait_for_pause_file(pause_path, check_interval):
  printed = False
  if pause_path and pause_path.is_file():
    if not printed:
      print('Execution paused..')
      printed = True
    time.sleep(check_interval)


class Thresholds:

  def __init__(self, args):
    self.values = {}
    for thres_name, thres_type in THRESHOLDS.items():
      thres_raw_value = getattr(args, thres_name)
      thres_value, thres_path = parse_file_or_value(thres_raw_value, thres_type)
      if thres_value is None and thres_path is None:
        self.values[thres_name] = None
      else:
        self.values[thres_name] = {'value':thres_value, 'path':thres_path, 'type':thres_type}

  def __getattr__(self, name):
    if self.values[name] is None:
      return None
    if self.values[name]['value'] is not None:
      return self.values[name]['value']
    elif self.values[name]['path'] is not None:
      return read_file(self.values[name]['path'], self.values[name]['type'])

  def subdivide_thres(self, main, specifics):
    """Use the value of threshold `main` as the default value for thresholds `specifics`."""
    for specific in specifics:
      if self.values[specific] is None:
        self.values[specific] = self.values[main].copy()


def parse_file_or_value(raw_value, coerce_type):
  if raw_value is None:
    return None, None
  try:
    return coerce_type(raw_value), None
  except ValueError:
    path = pathlib.Path(raw_value)
  if not path.is_file():
    fail(f'Error: Argument {raw_value!r} not a {coerce_type.__name__} or existing file.')
  return None, path


def get_idle(resource, min_node_size=None, chooser=max):
  assert resource in ('cpus', 'nodes'), resource
  assert hasattr(chooser, '__call__') or chooser in ('min', 'max'), chooser
  if not hasattr(chooser, '__call__'):
    if chooser == 'min':
      chooser = min
    elif chooser == 'max':
      chooser = max
  idle = 0
  best_cpus = None
  cmd = ('sinfo', '-h', '-p', 'general', '-t', 'idle,alloc', '-o', '%n %C')
  stdout = run_command(cmd, 'Error: Problem getting CPU usage info.')
  for line in stdout.splitlines():
    major_fields = line.split()
    if len(major_fields) != 2:
      continue
    node_name = major_fields[0]
    minor_fields = major_fields[1].split('/')
    if len(minor_fields) != 4:
      continue
    try:
      node_idle = int(minor_fields[1])
      node_size = int(minor_fields[3])
    except ValueError:
      continue
    if min_node_size is not None:
      if node_size < min_node_size:
        continue
    if resource == 'cpus':
      idle += node_idle
    elif resource == 'nodes':
      if node_idle == node_size:
        idle += 1
    elif resource == 'node':
      if best_cpus is None:
        best_cpus = node_idle
        idle = node_name
      else:
        result = chooser(node_idle, best_cpus)
        if result != best_cpus:
          best_cpus = node_idle
          idle = node_name
  return idle


def count_running_jobs(name=None):
  jobs = 0
  cmd = ('squeue', '-h', '-u', USER, '-t', 'running', '-o', '%j')
  stdout = run_command(cmd, 'Problem getting a list of running jobs.')
  for line in stdout.splitlines():
    if name is None:
      jobs += 1
    elif line == name:
      jobs += 1
  return jobs


def run_command(cmd, message):
  logging.info('Info: Running $ '+' '.join(map(str, cmd)))
  result = subprocess.run(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
  if result.returncode != 0:
    sys.stderr.write(str(result.stderr, 'utf8')+'\n')
    fail('Error: '+message)
  return str(result.stdout, 'utf8')


def read_file(path, coerce_type=None):
  with path.open('rt') as file:
    for line_raw in file:
      line = line_raw.strip()
      if coerce_type:
        return coerce_type(line)
      else:
        return line


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
