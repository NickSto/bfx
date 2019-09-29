#!/usr/bin/env python3
import argparse
import getpass
import logging
import pathlib
import subprocess
import sys
import time
assert sys.version_info.major >= 3, 'Python 3 required'

PARAMS = {
  'min_idle_nodes':int,
  'min_idle_cpus':int,
  'min_node_size':int,
  'min_node_size_cpus':int,
  'min_node_size_nodes':int,
  'max_jobs':int,
  'min_jobs':int,
  'cpus':int,
}
USER = getpass.getuser()
PAUSE_TIME = 10
USAGE = """$ %(prog)s [parameters]
       $ %(prog)s @path/to/args.txt"""
DESCRIPTION = """Determine whether we should keep launching slurm jobs, based on available
resources. Launch this and it will sleep until enough resources are available. Then, it will print
the name of a node with enough free CPUs to run your job."""
EPILOG = """You can also give arguments via a file. Just give the path to the file as an argument,
prepended by '@'. The file will be read, and each line will be interpreted as a regular argument to
this script. Remember to put literally every argument on its own line, so '--min-jobs 8' should be
given as two lines: '--min-jobs' and '8'."""


def make_argparser():
  parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG,
    fromfile_prefix_chars='@')
  parser.add_argument('-p', '--prefer', choices=('min', 'max'), default='min',
    help='Prefer nodes with either the most (max) or least (min) number of idle CPUs. '
      'Give --cpus to indicate how many CPUs the job requires. Only nodes with at least this '
      'many CPUs will be considered.')
  parser.add_argument('-C', '--cpus', type=int, default=1,
    help="How many CPUs are required by the job we're waiting to start. Default: %(default)s")
  parser.add_argument('-q', '--wait-for-job',
    help="Wait until the job with this name has begun. Useful if you just launched one and don't "
      "want to keep queueing jobs if they're not starting.")
  #TODO: Update help text with new, predictive algorithm.
  parser.add_argument('-n', '--min-idle-nodes', default=0,
    help='Keep this many nodes idle: wait if only this many + 1 are idle. Can give a number or a '
      'file containing the number (and nothing else).')
  parser.add_argument('-c', '--min-idle-cpus', default=0,
    help='Keep this many CPUs idle: wait if only this many + --cpus are idle. Can give a number or '
      'a file containing the number (and nothing else).')
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
    help='Always let yourself (try to) run at least this many jobs. Even if too few resources are '
      "available and you'd normally wait, keep going if fewer than this many jobs are running. "
      'In that case, this will exit but print nothing to stdout. Note: Does not override '
      '--pause-file. Can give a number or a file containing the number (and nothing else).')
  parser.add_argument('-P', '--pause-file', type=pathlib.Path,
    help='Wait if this file exists. Can be used as a manual pause button.')
  parser.add_argument('-i', '--check-interval', type=int, default=15,
    help='How many seconds to wait between checks for available resources. Default: %(default)s')
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

  params = Parameters(args)
  params.subdivide_param('min_node_size', ('min_node_size_cpus', 'min_node_size_nodes'))
  params.chooser = args.prefer

  if params.max_jobs is not None and params.min_jobs >= params.max_jobs:
    fail(f'Error: --min-jobs must be < --max-jobs ({params.min_jobs} >= {params.max_jobs}).')

  node = None
  wait = True
  last_reason = None
  while wait or node is None:
    wait = False
    paused = False
    reason_msg = None
    if args.wait_for_job and not count_running_jobs(args.wait_for_job):
      reason_msg = f'Waiting for job {job_name!r} to begin..'
      wait = True
    if count_running_jobs() >= params.max_jobs:
      reason_msg = f'Too many jobs running ({count_running_jobs()} >= {params.max_jobs})'
      wait = True
    if args.pause_file and args.pause_file.is_file():
      reason_msg = f'Execution paused: Pause file {args.pause_file} exists.'
      wait = True
      paused = True
    states = get_node_states()
    node = choose_node(
      states,
      args.cpus,
      min_idle_cpus=params.min_idle_cpus,
      min_idle_nodes=params.min_idle_nodes,
      min_node_cpus=params.min_node_cpus,
      chooser=params.prefer,
    )
    if count_running_jobs() < params.min_jobs and not paused:
      if wait or node is None:
        logging.info(
          f"You're running fewer than {params.min_jobs} jobs. Ignoring limits and continuing."
        )
      break
    if wait or node is None:
      if reason_msg and reason_msg != last_reason:
        logging.warning(reason_msg)
      last_reason = reason_msg
      time.sleep(args.check_interval)

  if node is not None:
    print(abbrev_node(node))


class Parameters:

  def __init__(self, args):
    self.values = {}
    for name, ptype in PARAMS.items():
      raw_value = getattr(args, name)
      value, path = parse_file_or_value(raw_value, ptype)
      if value is None and path is None:
        self.values[name] = None
      else:
        self.values[name] = {'value':value, 'path':path, 'type':ptype}

  def __getattr__(self, name):
    if self.values[name] is None:
      return None
    if self.values[name]['value'] is not None:
      return self.values[name]['value']
    elif self.values[name]['path'] is not None:
      return read_file(self.values[name]['path'], self.values[name]['type'])

  def subdivide_param(self, main, specifics):
    """Use the value of parameter `main` as the default value for parameters `specifics`."""
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


def get_node_states():
  states = {}
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
    states[node_name] = {'name':node_name, 'idle':node_idle, 'cpus':node_size}
  return states


def choose_node(
    states,
    job_cpus,
    min_idle_cpus=0,
    min_idle_nodes=0,
    min_node_size_cpus=1,
    min_node_size_idle=1,
    chooser=max
  ):
  """Choose a node to run the job on, if any.
  If the resources the job would consume would make them fall below the given thresholds, return
  `None`.
  `min_node_cpus`: Don't consider nodes with fewer than this many CPUs. Pretend they don't exist.
  `chooser`: Whether to prefer nodes with more or less free CPUs. `max` will make it prefer nodes
  with the most idle CPUs, spreading your jobs out across nodes. `min` will make it prefer nodes
  with fewer available CPUs (but still enough to run the job)."""
  chooser = get_chooser(chooser)
  idle_nodes, idle_cpus = count_idle_resources(
    states,
    min_node_size_cpus=min_node_size_cpus,
    min_node_size_idle=min_node_size_idle,
  )
  if idle_cpus - job_cpus < min_idle_cpus:
    return None
  if idle_nodes - 1 < min_idle_nodes:
    exclude_idle_nodes = True
  else:
    exclude_idle_nodes = False
  best_node = None
  for node in states.values():
    if node['cpus'] < min_node_size_cpus:
      continue
    if node['idle'] < job_cpus:
      continue
    if node['idle'] == node['cpus'] and exclude_idle_nodes:
      continue
    if best_node is None:
      best_node = node
    else:
      result = chooser(node['idle'], best_node['idle'])
      if result != best_node['idle']:
        best_node = node
  if best_node is None:
    return None
  else:
    return best_node['name']


def count_idle_resources(states, min_node_size_cpus=1, min_node_size_idle=1):
  idle_nodes = 0
  idle_cpus = 0
  for node in states.values():
    if node['cpus'] >= min_node_size_cpus:
      idle_cpus += node['idle']
    if node['cpus'] >= min_node_size_idle:
      if node['idle'] == node['cpus']:
        idle_nodes += 1
  return idle_nodes, idle_cpus


def get_chooser(chooser_raw):
  if hasattr(chooser_raw, '__call__'):
    return chooser_raw
  else:
    if chooser_raw == 'min':
      return min
    elif chooser_raw == 'max':
      return max
    else:
      fail(f'Error: Invalid chooser {chooser_raw!r}')


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


def abbrev_node(node_name):
  fields = node_name.split('.')
  return fields[0]


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
