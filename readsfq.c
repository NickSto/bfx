#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>

#define BUFFER_SIZE_DEFAULT 65535
const char *USAGE = "Usage: $ readsfq [options] reads.fq\n"
"       $ gunzip -c reads.fq.gz | readsfq [options]\n"
"This will count the number of reads in a FASTQ file, giving an accurate count\n"
"even for files with multi-line reads.\n"
"This will also report the range of quality score values, to help in determining\n"
"the format of the quality scores.\n"
"Options:\n"
"-q: Don't print anything to stderr.\n"
"-o [value]: What value to print to stdout. Default: \"reads\". Options:\n"
"   \"reads\":  The number of reads in the file.\n"
"   \"format\": The guessed format of the quality scores.\n"
"             Either \"sanger\" or \"solexa\".\n"
"-B [buffer_size]: Specify a file reading buffer size, in bytes. Default: 65535\n"
"                  WARNING: Lines longer than this will cause errors.\n"
"Note: The argument parsing is dumb. Don't combine options like \"-qo format\".";

//TODO: bool
typedef enum {
  HEADER, SEQ, PLUS, QUAL
} State;

typedef struct extremes {
  int max;
  int min;
} Extremes;

int line_is_empty(char *line);
void die(const char *message, ...);
int is_int(const char *int_str);
long count_chars(char *buffer, size_t buffer_size);
long count_chars_and_extremes(char *buffer, size_t buffer_size, Extremes *extremes);
char guess_quality_format(Extremes extremes, long num_reads);


int main(int argc, char *argv[]) {

  size_t buffer_size = BUFFER_SIZE_DEFAULT;

  // Read arguments
  FILE *infile = stdin;
  int quiet = 0;
  char output = 'r';
  char read_opt = '\0';
  int i;
  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      die(USAGE);
    } else if (strcmp(argv[i], "-q") == 0) {
      quiet = 1;
    } else if (strcmp(argv[i], "-o") == 0) {
      read_opt = 'o';
    } else if (strcmp(argv[i], "-B") == 0) {
      read_opt = 'B';
    } else if (read_opt == 'o') {
      if (strcmp(argv[i], "reads") == 0) {
        output = 'r';
      } else if (strcmp(argv[i], "format") == 0) {
        output = 'f';
      } else {
        die("Invalid -o output format \"%s\"", argv[i]);
      }
      read_opt = '\0';
    } else if (read_opt == 'B') {
      if (! is_int(argv[i])) {
        die("Invalid buffer size: \"%s\"", argv[i]);
      }
      buffer_size = atoi(argv[i]);
      read_opt = '\0';
    } else if (infile == stdin) {
      infile = fopen(argv[i], "r");
      if (errno) {
        die("\"%s\"", argv[i]);
      }
    } else {
      //TODO: allow any number of input files
      die("Can only process one file argument");
    }
  }

  int get_extremes = 1;
  if (quiet && output != 'f') {
    get_extremes = 0;
  }

  /*TODO: This assumes that there will be at least as many quality scores as there are sequence
   *      bases. According to Dan, we can't make that assumption.
   *      Then what do we do to tell when the quality lines have ended?
   *      Ideas for disambiguating:
   *      1. If len(qual) >= len(seq), it's a HEADER (If we've already seen enough
   *         quality values to cover the read, the QUAL lines must be over.)
   *      2. If the line plus the observed quality values so far is longer than the
   *         read, it must be a HEADER line.
   *      3. No FASTQ format uses space characters for quality scores, according to
   *         Wikipedia. If there's a space character in the line, say it's a HEADER line?
   *      But there could still conceivably be a read with truncated quality scores,
   *      followed by a HEADER line that contains no spaces and is short enough to not
   *      exceed the read length.
   *      Conclusion: Just check how BioPython does it:
   *      http://biopython.org/DIST/docs/api/Bio.SeqIO.QualityIO-pysrc.html
   *      Update: BioPython just throws an error if the number of bases != num of quality scores.
   */
  /* Notes on format requirements:
   * Empty lines are allowed, and ignored. If there's an empty line where a sequence or quality line
   * is expected, it's interpreted as 0 base/quality scores.
   * This means you can have zero-length reads.
   * Multi-line sequences (more than 4 lines per read) are allowed.
   * The number of quality scores must be >= the number of sequence bases. If there are missing
   * quality scores, this will likely fail with an error, but it could possibly succeed, giving
   * incorrect results.
   */

  char *line = malloc(buffer_size);

  Extremes extremes;
  extremes.max = 0;
  extremes.min = 256;
  long num_reads = 0;
  long seq_len = 0;
  long qual_len = 0;
  State state = HEADER;
  // fgets() reads a line at a time.
  char *result = fgets(line, buffer_size, infile);
  long line_num = 0;
  while (result != NULL) {
    line_num++;
    if (state == HEADER) {
      // Allow empty lines before the header.
      if (! line_is_empty(line)) {
        if (line[0] != '@') {
          die("Line %ld looked like a header line but does not start with \"@\".", line_num);
        }
        num_reads++;
        seq_len = 0;
        // Assume only 1 header line.
        state = SEQ;
      }
    } else if (state == SEQ) {
      if (line[0] == '+') {
        qual_len = 0;
        // End of sequence line comes when we see a line starting with "+".
        state = PLUS;
      } else {
        seq_len += count_chars(line, buffer_size);
      }
    } else if (state == PLUS || state == QUAL) {
      // If the state is PLUS, we already saw the "+" line on the last loop.
      // Assume there's only 1 "+" line, and assume we're now on a quality scores line.
      if (state == QUAL && line[0] == '@') {
        // If we're past the "first" quality scores line and we see one that starts with a "@",
        // that's very suspicious. Allow it, but raise a warning.
        fprintf(stderr, "Warning: Looking for more quality scores on line %ld but it starts with "
                        "\"@\".\nThis might be a header line and there were fewer quality scores "
                        "than bases.\n", line_num);
      }
      state = QUAL;
      if (get_extremes) {
        qual_len += count_chars_and_extremes(line, buffer_size, &extremes);
      } else {
        qual_len += count_chars(line, buffer_size);
      }
      if (qual_len >= seq_len) {
        // End of quality line comes once we've seen enough quality scores to match the sequence line.
        state = HEADER;
        if (qual_len > seq_len) {
          fprintf(stderr, "Warning on line %ld: Counted more quality scores than bases.\n", line_num);
        }
      }
    }
    result = fgets(line, buffer_size, infile);
  }

  char format_guess = '?';
  if (get_extremes) {
    format_guess = guess_quality_format(extremes, num_reads);
  }

  if (!quiet) {
    fprintf(stderr, "Quality score ascii range: %d (%c) to %d (%c)\n",
            extremes.min, (char)extremes.min, extremes.max, (char)extremes.max);
    switch (format_guess) {
      case 'S':
        fprintf(stderr, "Format: Very likely Sanger (offset 33).\n");
        break;
      case 'X':
        fprintf(stderr, "Format: Very likely Solexa (offset 64).\n");
        break;
      case 's':
        fprintf(stderr, "Format: Maybe Sanger? (offset 33)\n");
        break;
      case 'x':
        fprintf(stderr, "Format: Maybe Solexa? (offset 64)\n");
        break;
      case '?':
        fprintf(stderr, "Format: Unknown\n");
    }
  }

  if (output == 'r') {
    printf("%ld\n", num_reads);
  } else if (output == 'f') {
    switch (format_guess) {
      case 'S':
      case 's':
        printf("sanger\n");
        break;
      case 'X':
      case 'x':
        printf("solexa\n");
        break;
      default:
        printf("?\n");
    }
  }

  fclose(infile);
  return 0;
}


int line_is_empty(char *line) {
  switch (line[0]) {
    case '\0':
      // Empty string.
      return 1;
    case '\n':
      switch (line[1]) {
        case '\0':
          // Empty unix line.
          return 1;
      }
    case '\r':
      switch (line[1]) {
        case '\0':
          // Empty mac line.
          return 1;
        case '\n':
          if (line[2] == '\0') {
            // Empty dos line.
            return 1;
          }
      }
  }
  return 0;
}


// Count the number of content characters in the line.
// E.g. count the number of characters before any null, newline, or the end of the buffer.
long count_chars(char *buffer, size_t buffer_size) {
  long i = 0;
  while (i < buffer_size && buffer[i] != '\n' && buffer[i] != '\r' && buffer[i] != '\0') {
    i++;
  }
  return i;
}


long count_chars_and_extremes(char *buffer, size_t buffer_size, Extremes *extremes) {
  long i = 0;
  int value = 0;
  while (i < buffer_size && buffer[i] != '\n' && buffer[i] != '\r' && buffer[i] != '\0') {
    value = (int)buffer[i];
    if (value > extremes->max) {
      extremes->max = value;
    }
    if (value < extremes->min) {
      extremes->min = value;
    }
    i++;
  }
  return i;
}


// Estimate which quality score encoding it is based on the maximum and minimum ASCII value.
// Returns 'S' if it's confident it's Sanger, 's' if it's less confident,
// 'X' if it's confident it's Solexa, 'x' if it's less confident,
// and '?' if it's not sure enough to make a guess.
char guess_quality_format(Extremes extremes, long num_reads) {
  if (extremes.min < 59) {
    // Solexa values should never be below 59 (";", Solexa PHRED -5). If they are, it must be Sanger.
    return 'S';
  } else if (extremes.min >= 66) {
    if (extremes.max >= 76) {
      // If it never goes below 66 ("B", Solexa PHRED 2, Sanger PHRED 33),
      // and it ranges above 76 ("L", Solexa PHRED 12, Sanger PHRED 43),
      // then we can be pretty positive it's Solexa.
      if (num_reads > 100) {
        return 'X';
      } else {
        // But if we didn't see a ton of reads, we're less positive.
        return 'x';
      }
    } else {
      // If it stays between 66 and 76, that's a pretty weird and narrow range.
      // That's either PHRED 2 to 12 (Solexa) or 33 to 43 (Sanger).
      // Sanger might be barely more plausible, but not by much.
      return '?';
    }
  } else if (extremes.min >= 59) {
    if (extremes.max >= 76 || num_reads > 100) {
      // If the minimum is within the Solexa range and the maximum is above the reasonable Sanger
      // range, it's probably Solexa. Even if the maximum isn't out of range for Sanger, it's still
      // probably Solexa, since that's a really high minimum if it's Sanger.
      return 'x';
    } else {
      // However, if we just haven't seen many reads, then maybe the minimum is high just from a
      // small sample. Let's not say either way.
      // Note: It's pretty weird to have quality values all between 59-76.
      // That's either PHRED -5 to 12 (Solexa) or 26 to 43 (Sanger).
      return '?';
    }
  }
  return '?';
}


void die(const char *message, ...) {
  const char template[] = "Error: %s";
  size_t length = strlen(message) + strlen(template) - strlen("%s") + 1;
  // get argument list
  va_list data;
  va_start(data, message);
  // copy message into template
  char final_str[length];
  sprintf(final_str, template, message);
  // print to stderr and exit
  vfprintf(stderr, final_str, data);
  if (errno) {
    perror(" ");
  } else {
    fprintf(stderr, "\n");
  }
  exit(1);
}


/* Check if a string is a valid integer. Return 0 if invalid, 1 otherwise.
 * Checks if any character is outside the ASCII range '0' to '9', or if it's
 * greater than 10 digits (will overflow a 32-bit int).
 */
int is_int(const char *int_str) {
  size_t length = strlen(int_str);
  // If the string is empty or there are more than 10 digits (length of 2^32).
  // This won't catch all int overflows, but it won't false positive.
  if (length <= 0 || length > 10) {
    return 0;
  }
  int i;
  for (i = 0; i < length; i++) {
    // If the character's ASCII value is outside the "0" to "9" range.
    if (int_str[i] < 48 || int_str[i] > 57) {
      return 0;
    }
  }
  return 1;
}
