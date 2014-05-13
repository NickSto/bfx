#!/usr/bin/awk
# original author: Nick Stoler
# Print the number of instances of single reads, read pairs, read triples, etc.
# Input should be a SAM with no header.
# N.B.: INPUT MUST BE SORTED BY READ NAME

BEGIN {
  OFS = "\t"
  FS = "\t"
  occurrences = 1
}

last == $1 {
  occurrences++
}
last != $1 && NR > 1 {
  totals[occurrences]++
  occurrences = 1
}
{
  last = $1
}

END {
  totals[occurrences]++
  for (occurrences in totals) {
    printf "%d\t%d\n", occurrences, totals[occurrences]
  }
}