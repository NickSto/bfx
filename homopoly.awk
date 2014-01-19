#!/usr/bin/env bioawk
# Print the location, length, and base for all homopolymers above a given
# length in a FASTA sequence.
# USAGE: $ bioawk -f homopoly.awk -c fastx -v thres=4 ~/data/chrM-rCRS.fa

BEGIN {
  MAX_LEN = 5000000; # 5Mb = ~500Mb of memory
  if (! thres) {
    thres = 3;
  }
}

{
  len = length($seq);
  if (len > MAX_LEN) {
    print $name, len, "Sequence too long. Skipping.";
  } else {
    split($seq, bases, "");
    for (i = 1; i <= len; i++) {
      base = toupper(bases[i]);
      if (base == last) {
        tot++;
      } else {
        if (tot >= thres) {
          print start, tot, last;
        }
        tot = 1;
        start = i;
      }
      last = base;
    }
  }
}
