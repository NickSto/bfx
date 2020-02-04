#!/usr/bin/env python3
import argparse
import configparser
import getpass
import logging
import pathlib
import subprocess
import sys
import time
assert sys.version_info.major >= 3, 'Python 3 required'

def boolish(raw):
  if raw.lower() in ('true', '1'):
    return True
  elif raw.lower() in ('false', '0'):
    return False
  else:
    return None

PARAMS = {
  'min_idle_nodes': {'type':int, 'default':0},
  'min_idle_cpus': {'type':int, 'default':0},
  'min_node_size': {'type':int, 'default':1},
  'min_node_size_cpus': {'type':int, 'fallback':'min_node_size'},
  'min_node_size_nodes': {'type':int, 'fallback':'min_node_size'},
  'max_jobs': {'type':int},
  'min_jobs': {'type':int, 'default':0},
  'prefer': {'type':str, 'default':'min'},
  'cpus': {'type':int, 'default':1},
  'mem': {'type':int, 'default':0},
  'stop': {'type':boolish, 'default':False},
}
PARAM_TYPES = {name:meta['type'] for name, meta in PARAMS.items()}
UNITS = {'B':1, 'K':1024, 'M':1024**2, 'G':1024**3, 'T':1024**4}
USER = getpass.getuser()
PAUSE_TIME = 10
USAGE = """$ %(prog)s [parameters]
       $ %(prog)s -c config.ini"""
DESCRIPTION = """Determine whether we should keep launching slurm jobs, based on available
resources. Launch this and it will sleep until enough resources are available. Then, it will print
the name of a node with enough free CPUs to run your job."""


def make_argparser():
  parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION, add_help=False)
  options = parser.add_argument_group('Options')
  options.add_argument('-c', '--config', type=pathlib.Path,
    help='A config file to set all the parameters below. This will be read after every pause, so '
      'you can update it while this script is running and it will change its behavior.')
  options.add_argument('-q', '--wait-for-job',
    help="Wait until the job with this name has begun. Useful if you just launched one and don't "
      "want to keep queueing jobs if they're not starting.")
  options.add_argument('-Q', '--wait-for-job-prefix',
    help='Same as --wait-for-job, but accept any job whose name starts with this string.')
  options.add_argument('-P', '--pause-file', type=pathlib.Path,
    help='Wait if this file exists. Can be used as a manual pause button.')
  options.add_argument('-i', '--check-interval', type=int, default=15,
    help='How many seconds to wait between checks for available resources. Default: %(default)s')
  options.add_argument('--mock-sinfo', type=pathlib.Path)
  options.add_argument('-h', '--help', action='help',
    help='Print this argument help text and exit.')
  params = parser.add_argument_group('Parameters')
  params.add_argument('-C', '--cpus', type=int,
    help="How many CPUs are required by the job we're waiting to start. Default: "
      +str(PARAMS['cpus']['default']))
  params.add_argument('-M', '--mem', type=bytes_spec,
    help="How much memory is required by the job we're waiting to start. You can give a number "
      'ending "B", "K", "M", "G", or "T", case-insensitive. Default: '
      +str(PARAMS['mem']['default']))
  #TODO: Update help text with new, predictive algorithm.
  params.add_argument('-n', '--min-idle-nodes',
    help='Keep this many nodes idle: wait if only this many + 1 are idle.')
  params.add_argument('-u', '--min-idle-cpus',
    help='Keep this many CPUs idle: wait if only this many + --cpus are idle.')
  params.add_argument('-s', '--min-node-size',
    help="Minimum node size when counting available resources for the above thresholds. Don't "
      'consider nodes with fewer than this many CPUs.')
  params.add_argument('--min-node-size-cpus',
    help='Same as --min-node-size, but only for when counting idle cpus for --min-idle-cpus')
  params.add_argument('--min-node-size-nodes',
    help='Same as --min-node-size, but only for when counting idle nodes for --min-idle-nodes')
  params.add_argument('-J', '--max-jobs',
    help="Don't let yourself have more than this many jobs running at once. If you have this many "
      'jobs running, wait.')
  params.add_argument('-j', '--min-jobs',
    help='Always let yourself (try to) run at least this many jobs. Even if too few resources are '
      "available and you'd normally wait, keep going if fewer than this many jobs are running. "
      'In that case, this will exit but print nothing to stdout. Note: Does not override '
      '--pause-file.')
  params.add_argument('-p', '--prefer', choices=('min', 'max'),
    help='Prefer nodes with either the most (max) or least (min) number of idle CPUs. '
      'Give --cpus to indicate how many CPUs the job requires. Only nodes with at least this '
      'many CPUs will be considered. Default: '+str(PARAMS['prefer']['default']))
  log = parser.add_argument_group('Logging')
  log.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = log.add_mutually_exclusive_group()
  volume.add_argument('--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  if args.config and not args.config.is_file():
    logging.warning(f'Warning: --config file {str(args.config)!r} not found!')

  if args.wait_for_job:
    if args.wait_for_job_prefix:
      fail('Error: Cannot give both --wait-for-job and --wait-for-job-prefix.')
    wait_for = args.wait_for_job
    prefixed = False
  else:
    wait_for = args.wait_for_job_prefix
    prefixed = True

  params = Parameters(args=args, config=args.config)

  if params.max_jobs is not None and params.min_jobs >= params.max_jobs:
    fail(f'Error: --min-jobs must be < --max-jobs ({params.min_jobs} >= {params.max_jobs}).')

  node = None
  wait = True
  last_reason = None
  while wait or node is None:
    wait = False
    paused = False
    reason_msg = None
    if wait_for and not count_running_jobs(name=wait_for, prefixed=prefixed):
      reason_msg = f'Waiting for job {wait_for!r} to begin..'
      wait = True
    if params.max_jobs and count_running_jobs() >= params.max_jobs:
      reason_msg = f'Too many jobs running ({count_running_jobs()} >= {params.max_jobs})'
      wait = True
    if args.pause_file and args.pause_file.is_file():
      reason_msg = f'Execution paused: Pause file {args.pause_file} exists.'
      wait = True
      paused = True
    if args.mock_sinfo:
      states = read_mock_sinfo(args.mock_sinfo)
    else:
      states = get_node_states()
    node = choose_node(
      states,
      params.cpus,
      params.mem,
      min_idle_cpus=params.min_idle_cpus,
      min_idle_nodes=params.min_idle_nodes,
      min_node_size_cpus=params.min_node_size_cpus,
      min_node_size_nodes=params.min_node_size_nodes,
      chooser=params.prefer,
    )
    if node is None and reason_msg is None:
      reason_msg = 'No node currently fits the given constraints..'
    if params.min_jobs and count_running_jobs() < params.min_jobs and not paused:
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
      if args.config:
        params.update_with_config(args.config)
    if params.stop:
      logging.warning(f'Instructed to stop by {args.config}.')
      node = 'STOP'
      break

  if node is not None:
    print(abbrev_node(node))


class Parameters:

  def __init__(self, args=None, config=None):
    # Initialize params.
    self.values = {}
    for name, meta in PARAMS.items():
      self.values[name] = None
    # Update with config file, if any.
    if config:
      self.update_with_config(config)
    # Update with args, if any.
    if args:
      self.update_with_args(args)
    self.set_defaults()

  def update_with_config(self, config_path):
    config = read_config_section(config_path, 'params', types=PARAM_TYPES)
    for name, value in config.items():
      if value is not None:
        self.values[name] = value

  def update_with_args(self, args):
    for name, meta in PARAMS.items():
      raw_value = getattr(args, name, None)
      if raw_value is not None:
        value = meta['type'](raw_value)
        self.values[name] = value

  def __getattr__(self, name):
    value = self.values.get(name)
    if value is not None:
      return value
    else:
      return None

  def set_defaults(self):
    """Set any unset params to their defaults (or the defaults of their fallbacks)."""
    for name, meta in PARAMS.items():
      if self.values[name] is not None:
        continue
      if 'default' in meta:
        self.values[name] = meta['default']
      elif 'fallback' in meta:
        fallback = meta['fallback']
        self.values[name] = self.values[fallback]

  def subdivide_param(self, main, specifics):
    """Use the value of parameter `main` as the default value for parameters `specifics`.
    Deprecated."""
    for specific in specifics:
      if self.values[specific] is None:
        self.values[specific] = self.values[main]


def read_config_section(config_path, section, types=None):
  data = {}
  config = configparser.ConfigParser(interpolation=None)
  try:
    config.read(config_path)
    for key, raw_value in config.items(section):
      if types and key in types:
        value = types[key](raw_value)
      else:
        value = raw_value
      data[key] = value
  except configparser.Error as error:
    fail(f'Invalid config file format in {config_path!r}: {error}')
  return data


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
  cmd = ('sinfo', '-h', '-p', 'general', '-t', 'idle,alloc', '-o', '%n %C %e')
  stdout = run_command(cmd, 'Error: Problem getting CPU usage info.')
  for line in stdout.splitlines():
    major_fields = line.split()
    if len(major_fields) != 3:
      logging.warning(
        f'Warning: sinfo line has wrong number of fields ({len(major_fields)}): {line!r}'
      )
      continue
    node_name = major_fields[0]
    minor_fields = major_fields[1].split('/')
    mem_str = major_fields[2]
    if len(minor_fields) != 4:
      logging.warning(
        f'Warning: sinfo line has wrong number of cpu fields ({len(minor_fields)}): {line!r}'
      )
      continue
    try:
      node_idle = int(minor_fields[1])
      node_size = int(minor_fields[3])
    except ValueError:
      logging.warning(
        f'Warning: sinfo line has invalid cpu fields ({minor_fields[1]!r} or {minor_fields[3]!r}): '
        f'{line!r}'
      )
      continue
    try:
      mem = int(mem_str) * 1024**2
    except ValueError:
      logging.warning(f'Warning: sinfo line has invalid mem field ({mem_str!r}): {line!r}')
      continue
    states[node_name] = {'name':node_name, 'idle':node_idle, 'cpus':node_size, 'mem':mem}
  return states


def choose_node(
    states,
    job_cpus,
    job_mem,
    min_idle_cpus=0,
    min_idle_nodes=0,
    min_node_size_cpus=1,
    min_node_size_nodes=1,
    chooser=max
  ):
  """Choose a node to run the job on, if any.
  If the resources the job would consume would make them fall below the given thresholds, return
  `None`.
  `chooser`: Whether to prefer nodes with more or less free CPUs. `max` will make it prefer nodes
  with the most idle CPUs, spreading your jobs out across nodes. `min` will make it prefer nodes
  with fewer available CPUs (but still enough to run the job)."""
  chooser = get_chooser(chooser)
  idle_nodes, idle_cpus = count_idle_resources(
    states,
    min_node_size_cpus=min_node_size_cpus,
    min_node_size_nodes=min_node_size_nodes,
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
    if node['mem'] < job_mem:
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


def count_idle_resources(states, min_node_size_cpus=1, min_node_size_nodes=1):
  idle_nodes = 0
  idle_cpus = 0
  for node in states.values():
    if node['cpus'] >= min_node_size_cpus:
      idle_cpus += node['idle']
    if node['cpus'] >= min_node_size_nodes:
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


def count_running_jobs(name=None, prefixed=False):
  jobs = 0
  cmd = ('squeue', '-h', '-u', USER, '-t', 'running', '-o', '%j')
  stdout = run_command(cmd, 'Problem getting a list of running jobs.')
  for line in stdout.splitlines():
    if name is None:
      jobs += 1
    else:
      if prefixed:
        if line.startswith(name):
          jobs += 1
      else:
        if line == name:
          jobs += 1
  return jobs


def abbrev_node(node_name):
  fields = node_name.split('.')
  return fields[0]


def bytes_spec(bytes_str):
  quantity_str = bytes_str[:len(bytes_str)-1]
  unit = bytes_str[-1].upper()
  # Possible ValueError to be caught by caller.
  quantity = int(quantity_str)
  if unit not in UNITS:
    raise ValueError(f'Invalid unit in byte amount {bytes_str!r}.')
  multiplier = UNITS[unit]
  return quantity * multiplier


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


def read_mock_sinfo(path):
  """Format: 4 space-delimited columns: node, total CPUs, free CPUs, free mem (GB)."""
  states = {}
  with path.open() as sinfo_file:
    for line in sinfo_file:
      fields = line.split()
      name = fields[0]
      cpus = int(fields[1])
      idle = int(fields[2])
      mem = int(fields[3]) * 1024**3
      states[name] = {'name':name, 'cpus':cpus, 'idle':idle, 'mem':mem}
  return states


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except (BrokenPipeError, KeyboardInterrupt):
    pass
