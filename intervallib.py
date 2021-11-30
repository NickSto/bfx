"""Interval operations."""

def read_intervals(intervals_file):
  """Read an intervals file.
  Format: Tab-delimited, 1st column is start coordinate, 2nd column is end coordinate (1-based)."""
  for line_raw in intervals_file:
    fields = line_raw.rstrip('\r\n').split('\t')
    start_str, end_str = fields[:2]
    start, end = int(start_str), int(end_str)
    yield (start, end)


def simplify_intervals(intervals):
  """Combine overlapping intervals, returning a list of regions covered by at least one interval.
  `intervals` must be in the same format as for `find_interval()`, except the intervals do not have
    to be sorted and they can be overlapping."""
  coverage = []
  last_interval = None
  for interval in sorted(intervals, key=lambda interval: interval[0]):
    if last_interval is None:
      coverage.append(interval)
    else:
      last_start, last_end = last_interval # pylint: disable=unpacking-non-sequence
      start, end = interval
      if last_end < start:
        coverage.append(interval)
      else:
        # This interval overlaps the previous one. Combine them and replace the previous interval
        # with the combination.
        coverage.pop()
        coverage.append((last_start, end))
    last_interval = interval
  return coverage


def find_interval(coord, intervals):
  """Find the first interval that contains the coordinate `coord`.
  This uses a binary search, so the time taken is O(log(n)).
  Arguments:
  `coord`: The coordinate to search with.
  `intervals`: A sequence of intervals, each of which is a sequence of two integers:
    the start and end coordinates. The sequence must not contain any intervals which overlap,
    and it must be sorted by start coordinate. In other words: the start of each interval must
    be *after* the end of the previous one.
  Returns:
  The interval which contains `coord` or `None` if no match is found. The interval start and end
  are considered inclusive, so if `start == coord` or `coord == end`, the interval is a match."""
  if not intervals:
    return None
  i = find_interval_i(coord, intervals)
  if i == -1:
    return None
  else:
    interval = intervals[i]
    if interval[0] <= coord <= interval[1]:
      return interval
    else:
      return None


def find_intervals(start, end, intervals):
  """Find all intervals between `start` and `end`.
  `intervals` must be in the same format as for `find_interval()`.
  This returns a list of intervals which overlap the interval `(start, end)`. This includes any
  interval which contains `start` OR `end`."""
  results = []
  if not intervals:
    return results
  i = find_interval_i(start, intervals)
  if i == -1:
    i = 0
  interval = intervals[i]
  started = False
  while interval[0] <= end:
    if not started and interval[1] >= start:
      started = True
    if started:
      results.append(interval)
    i += 1
    if i >= len(intervals):
      break
    interval = intervals[i]
  return results


def find_interval_i(coord, intervals, subset_start=0, subset_end=None):
  """Find the index of the interval containing the coordinate, or is closest to the point where such
  an interval would be.
  Arguments:
  `coord`: Same as in `find_interval()`.
  `intervals`: Same as in `find_interval()`.
  Returns:
  An index into the `intervals` sequence. If there is an interval containing `coord`, this will be
  the index of that interval. If no interval contains `coord`, this will return the index of the
  interval just to the "left" of `coord`. If `coord` is before the first interval, this will return
  -1. WARNING: If `intervals` is empty (there are no intervals), this returns `0`, which will raise
  an IndexError if used as an index for `intervals`."""
  if subset_end is None:
    subset_end = len(intervals) - 1
  num_intervals = subset_end - subset_start + 1
  if num_intervals <= 1:
    # Only one (or fewer) interval left.
    if subset_start == subset_end == 0 and coord < intervals[0][0]:
      # Special case: `coord` is before all intervals.
      return -1
    else:
      return subset_start
  i = subset_start + (num_intervals // 2)
  start, end = intervals[i]
  if start <= coord:
    if coord <= end:
      # We have a match!
      return i
    else:
      # `coord` is after this interval. Search the right half.
      return find_interval_i(coord, intervals, subset_start=i, subset_end=subset_end)
  else:
    # `coord` is before this interval. Search the left half.
    return find_interval_i(coord, intervals, subset_start=subset_start, subset_end=i-1)


def intersect_intervals(query_interval, intervals):
  """Get the intersection of the query interval with a list of other intervals.
  `intervals` must be in the same format as for `find_interval()`.
  Returns a new list of intervals: the portions of `intervals` which are inside `query_interval`."""
  start, end = query_interval
  results = []
  for interval in find_intervals(start, end, intervals):
    intersect_start = max(interval[0], start)
    intersect_end = min(interval[1], end)
    intersect = (intersect_start, intersect_end)
    assert intersect_start <= intersect_end, (start, end, interval, intersect)
    results.append(intersect)
  return results
