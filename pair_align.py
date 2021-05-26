#!/usr/bin/env python3
import argparse
import logging
import pathlib
import sys
try:
  from Bio.Align import PairwiseAligner
except ImportError:
  print(
    'Error importing Bio.Align.PairwiseAligner. Make sure BioPython is properly installed.',
    file=sys.stderr
  )
  raise
try:
  from Bio.Align import substitution_matrices
except ImportError:
  substitution_matrices = None
  print(
    'Error loading `Bio.Align.substitution_matrices`. Using shim instead. Install BioPython '
    '1.75 or later to use the native version.', file=sys.stderr
  )

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
DEFAULT_MATRIX = SCRIPT_DIR/'sub-matrix.txt'
DEFAULT_GAP_OPEN = 0
DEFAULT_GAP_EXTEND = 2.5
DESCRIPTION = """Align two sequences with BioPython's PairwiseAligner."""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  options.add_argument('seq1')
  options.add_argument('seq2')
  options.add_argument('-m', '--matrix', type=pathlib.Path,
    help=f'Path to a substitution matrix file. Default: {DEFAULT_MATRIX}')
  options.add_argument('-g', '--global', dest='scope', default='local', action='store_const',
    const='global',
    help='Do a global alignment, instead of a local one (the default).')
  options.add_argument('-o', '--gap-open', type=float, default=DEFAULT_GAP_OPEN,
    help='Gap open penalty. This will be subtracted from the alignment score at the start of each '
      'gap. Default: %(default)s')
  options.add_argument('-e', '--gap-extend', type=float, default=DEFAULT_GAP_EXTEND,
    help='Gap extension penalty. This will be subtracted from the alignment score for every base '
      'of the gap. Default: %(default)s')
  options.add_argument('-t', '--trim', action='store_true',
    help='Trim non-matching ends from the resulting alignment.')
  options.add_argument('-h', '--help', action='help',
    help='Print this argument help text and exit.')
  logs = parser.add_argument_group('Logging')
  logs.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = logs.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  if args.matrix:
    matrix_path = args.matrix
  else:
    matrix_path = DEFAULT_MATRIX

  alignment = align(
    args.seq1, args.seq2, scope=args.scope, gap_open=args.gap_open, gap_extend=args.gap_extend,
    matrix_file=matrix_path, trim=args.trim,
  )
  if alignment is None:
    logging.warning('Warning: No satisfactory alignment found.')
  else:
    print(alignment)


def align(seq1, seq2, **kwargs):
  aligner = Aligner(**kwargs)
  return aligner.align(seq1, seq2)


class Aligner:
  def __init__(
    self, scope='local', gap_open=DEFAULT_GAP_OPEN, gap_extend=DEFAULT_GAP_EXTEND,
    matrix_file=DEFAULT_MATRIX, matrix=None, trim=False,
  ):
    self._aligner = PairwiseAligner()
    if matrix is None:
      if not matrix_file:
        raise RuntimeError(f'Must give either a `matrix` or `matrix_file` to {type(self).__name__}')
      if substitution_matrices is None:
        matrix = substitution_matrix_loader(matrix_file)
      else:
        matrix = substitution_matrices.read(matrix_file)
    self.trim = trim
    self.scope = scope
    self._gap_extend = gap_extend
    self._gap_open = gap_open
    self.gap_extend = gap_extend
    self.gap_open = gap_open
    self._aligner.substitution_matrix = matrix

  @property
  def algorithm(self):
    return self._aligner.algorithm

  @property
  def scope(self):
    return self._aligner.mode
  @scope.setter
  def scope(self, value):
    self._aligner.mode = value

  @property
  def matrix(self):
    return self._aligner.substitution_matrix

  # The gap penalty values are actually *added* to the alignment score, so if you're using the
  # intuitive sense of the term "penalty" (a value to be subtracted from the score), you have to
  # negate it.
  # So  the PairwiseAligner.open_gap_score = -(gap_open + gap_extend)
  # and the PairwiseAligner.extend_gap_score = -gap_extend
  # Also, BioPython takes `open_gap_score` to mean the entirety of the penalty applied to the
  # first base of the gap. The `extend_gap_score` will be applied to subsequent bases, but not the
  # first base. So if you're using the more common definitions of "gap open" and "gap extension"
  # penalties, the `open_gap_score` must be the sum of your gap open and gap extension penalties.
  # See https://biopython.readthedocs.io/en/latest/chapter_align.html#affine-gap-scores for where
  # they define how they use the values.
  @property
  def gap_extend(self):
    return self._gap_extend
  @gap_extend.setter
  def gap_extend(self, value):
    self._gap_extend = value
    self._aligner.extend_gap_score = -self._gap_extend
    self._aligner.open_gap_score = -(self._gap_open + self._gap_extend)

  @property
  def gap_open(self):
    return self._gap_open
  @gap_open.setter
  def gap_open(self, value):
    self._gap_open = value
    self._aligner.open_gap_score = -(self._gap_open + self._gap_extend)

  def align(self, seq1, seq2):
    bp_alignments = self._aligner.align(seq1, seq2)
    if len(bp_alignments) > 0:
      alignment = Alignment(bp_alignment=bp_alignments[0], aligner=self)
      if self.trim:
        return alignment.trimmed()
      else:
        return alignment
    else:
      return None

  def aligns(self, seq1, seq2):
    bp_alignments = self._aligner.align(seq1, seq2)
    alignments = (Alignment(bp_alignment=aln, aligner=self) for aln in bp_alignments)
    if self.trim:
      return (aln.trimmed() for aln in alignments)
    else:
      return alignments


class Alignment:
  def __init__(
    self, target=None, matches_str=None, query=None, score=None, bp_alignment=None, aligner=None,
  ):
    if bp_alignment is None:
      if target is None or matches_str is None or query is None:
        raise ValueError(
          "Must either specify every attribute's value or give a BioPython PairwiseAlignment object."
        )
      self.target, self.matches_str, self.query, self.score = target, matches_str, query, score
    else:
      if not (target is None and matches_str is None and query is None and score is None):
        raise ValueError(
          "Must not specify individual attribute values when giving a BioPython PairwiseAlignment "
          "object."
        )
      self.target, self.matches_str, self.query = str(bp_alignment).splitlines()
      self.score = bp_alignment.score
    self.aligner = aligner

  def trimmed(self, thres=0.5, thres_score=None):
    """Get a "trimmed" version of the alignment where non-matching tails are removed.
    This is different from a local alignment and can be applied to a global one. It is a post-
    processing step that is applied after an alignment (even a global one).
    A "tail" is defined as sequence before the first "good match" in the alignment or after the last
    "good match". A "good match" is a non-gap position with a match score at or above the threshold:
    If `thres_score` is given, this is the minimum score needed. Any position with a match score
    (from the substitution matrix) lower than this is considered a bad match.
    Otherwise, the `thres_score` is derived from `thres`. This is a fraction between 0 and 1 which
    expresses the `thres_score` in terms of the range between the minimum and maximum scores in the
    substitution matrix. Specifically, `thres_score = thres * (max_score-min_score) + min_score`.
    """
    if thres_score is None:
      min_score = min(self.aligner.matrix.values())
      max_score = max(self.aligner.matrix.values())
      score_range = max_score - min_score
      thres_score = thres * score_range + min_score
    start = end = None
    for i, (base1, base2) in enumerate(zip(self.target, self.query)):
      # Is it a good enough match?
      good_match = True
      try:
        score = self.aligner.matrix[(base1,base2)]
      except IndexError:
        good_match = False
      else:
        if score < thres_score:
          good_match = False
      # If we're in good match territory, update the edge indices.
      if good_match:
        if start is None:
          start = i
        end = i
    if start is None or end is None:
      return None
    target = self.target[start:end+1]
    matches_str = self.matches_str[start:end+1]
    query = self.query[start:end+1]
    return Alignment(target=target, matches_str=matches_str, query=query, aligner=self.aligner)

  def __str__(self):
    """Print a human-readable representation of the alignment."""
    return '\n'.join((self.target, self.matches_str, self.query))


def substitution_matrix_loader(matrix_path):
  x_bases = None
  matrix = {}
  with matrix_path.open() as matrix_file:
    for line_num, line_raw in enumerate(matrix_file,1):
      line = line_raw.strip()
      if line.startswith('#'):
        continue
      if x_bases is None:
        x_bases = line.split()
        continue
      y_base, *score_strs = line.split()
      scores = [float(s) for s in score_strs]
      if len(scores) != len(x_bases):
        raise ValueError(
          f'Line {line_num} has {len(scores)} scores but the header has {len(x_bases)} bases.'
        )
      for x_base, score in zip(x_bases, scores):
        matrix[(x_base,y_base)] = score
  return matrix


def fail(message):
  logging.critical(f'Error: {message}')
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception(message)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
