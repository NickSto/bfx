#!/usr/bin/env python3
import argparse
import difflib
import logging
import pathlib
import subprocess
import sys
import types
assert sys.version_info.major >= 3, 'Python 3 required'

TESTS_DIR = pathlib.Path(__file__).resolve().parent
ROOT_DIR = TESTS_DIR.parent
DESCRIPTION = """"""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('tests', metavar='test_name', nargs='*',
    help='The tests to run.')
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

  # Create the dicts holding all defined tests.

  meta_tests = get_objects_diff(GlobalsInitial, GlobalsAfterMeta)
  meta_tests['all'] = lambda name: run_test_group(simple_tests)
  meta_tests['active'] = lambda name: run_test_group(active_tests)
  meta_tests['inactive'] = lambda name: run_test_group(inactive_tests)

  active_tests = get_objects_diff(GlobalsAfterMeta, GlobalsAfterActive)
  inactive_tests = get_objects_diff(GlobalsAfterActive, GlobalsAfterInactive)
  simple_tests = add_dicts(active_tests, inactive_tests)

  all_tests = add_dicts(meta_tests, simple_tests)

  if not args.tests:
    print('Meta tests:')
    for test_name in meta_tests.keys():
      print('  '+test_name)
    print('Active tests:')
    for test_name in active_tests.keys():
      print('  '+test_name)
    print('Inactive tests:')
    for test_name in inactive_tests.keys():
      print('  '+test_name)
    return 1

  unknown_tests = []
  for test_name in args.tests:
    if test_name not in all_tests:
      unknown_tests.append(test_name)
  if unknown_tests:
    fail('Error: Test(s) "{}" unrecognized.'.format('", "'.join(unknown_tests)))

  for test_name in args.tests:
    test_fxn = all_tests[test_name]
    test_fxn(test_name)

GlobalsInitial = globals().copy()


##### Meta tests #####

GlobalsAfterMeta = globals().copy()


##### Active tests #####

def getreads_smoke(test_name):
  script_name = 'getreads.py'
  script = ROOT_DIR / script_name
  test_pairs = (
    ('smoke.fq',  'smoke.fq.out'),
    ('smoke.fa',  'smoke.fa.out'),
    ('smoke.txt', 'smoke.txt.out'),
    ('smoke.tsv', 'smoke.tsv.out'),
    ('smoke.sam', 'smoke.sam.out'),
  )
  for input_name, output_name in test_pairs:
    print('{} ::: {} ::: {}\t'.format(test_name, script_name, input_name), end='')
    result = run_command((script, TESTS_DIR/input_name))
    expected = read_file(TESTS_DIR/output_name)
    if result != expected:
      print('FAILED')
      for line in trimmed_diff(expected.splitlines(), result.splitlines()):
        print(line)
    else:
      print('success')


GlobalsAfterActive = globals().copy()


##### Inactive tests #####

GlobalsAfterInactive = globals().copy()


##### Helper functions #####


def add_dicts(*dicts):
  combined = {}
  for d in dicts:
    for key, value in d.items():
      combined[key] = value
  return combined


def run_test_group(test_group):
  for test_name, test_fxn in test_group.items():
    test_fxn(test_name)


def get_objects_diff(objects_before, objects_after, object_type=types.FunctionType):
  diff = {}
  for name, obj in objects_after.items():
    if name not in objects_before and isinstance(obj, object_type):
      diff[name] = obj
  return diff


def run_command(command):
  try:
    output = subprocess.check_output(command, stderr=subprocess.DEVNULL)
  except OSError as error:
    logging.error('Error: {}'.format(error))
    return None
  except subprocess.CalledProcessError as error:
    logging.error('Error: {}'.format(error))
    return None
  return str(output, 'utf8')


def read_file(path):
  try:
    with path.open('r') as file:
      return file.read()
  except OSError:
    logging.error('Error: {}'.format(error))
    return None


def trimmed_diff(lines1, lines2, lineterm=''):
  """Get a trimmed diff.
  Input lines should be newline-free."""
  diff_lines = difflib.unified_diff(lines1, lines2, n=1, fromfile='a', tofile='b',
                                    fromfiledate='c', tofiledate='d', lineterm=lineterm)
  header_line = 0
  for line in diff_lines:
    if header_line == 0 and line == '--- a\tc'+lineterm:
      header_line = 1
    elif header_line == 1 and line == '+++ b\td'+lineterm:
      header_line = 2
    elif header_line == 2:
      header_line = None
    if header_line is None:
      yield line


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
