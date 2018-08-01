#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>

/* DISCLAIMER:
 * This is a My First C Program, and naiveley coded, but hopefully relatively
 * bug-free. I was just happy it worked.
*/

#define BUFFER_SIZE_DEFAULT 65535
const char *USAGE = "Usage: $ readsfq [options] reads.fq\n"
"       $ gunzip -c reads.fq.gz | readsfq [options]\n"
"This will count the number of reads in a FASTQ file, giving an accurate count\n"
"even for files with multi-line reads.\n"
"Options:\n"
"-B [buffer size]: Specify a file reading buffer size, in bytes. Default: 65535";

void die(const char *message, ...);
int is_int(const char *int_str);

//TODO: bool
typedef enum {
  FIRST, NAME, SEQ, PLUS, QUAL
} State;

int main(int argc, char *argv[]) {

  size_t buffer_size = BUFFER_SIZE_DEFAULT;

  // Read arguments
  FILE *infile = stdin;
  char read_opt = '\0';
  int i;
  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      die(USAGE);
    } else if (strcmp(argv[i], "-B") == 0) {
      read_opt = 'B';
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

  char *buffer = malloc(buffer_size);

  int num_reads = 0;
  State line_type = FIRST;
  short first_char = 1;
  //TODO: fgets() reads a line at a time.
  size_t bytes_read = fread(buffer, 1, buffer_size, infile);
  while (bytes_read) {
    char chr;
    int i;
    for (i = 0; i < bytes_read; i++) {
      chr = buffer[i];
      if (chr == '\n' || chr == '\r') {
        first_char = 1;
      } else if (first_char) {
        switch (chr) {
          case '@':
            if (line_type == FIRST || line_type == QUAL) {
              /*TODO: Pretty sure this is wrong. If the line after a QUAL line starts with a '@',
               *      it could be either another QUAL line or a NAME line. According to Dan, we
               *      can't assume the quality scores are as long as the sequence, since the quality
               *      scores can be truncated.
               *      Ideas for disambiguating:
               *      1. If len(qual) >= len(seq), it's a NAME (If we've already seen enough
               *         quality values to cover the read, the QUAL lines must be over.)
               *      2. If the line plus the observed quality values so far is longer than the
               *         read, it must be a NAME line.
               *      3. No FASTQ format uses space characters for quality scores, according to
               *         Wikipedia. If there's a space character in the line, say it's a NAME line?
               *      But there could still conceivably be a read with truncated quality scores,
               *      followed by a NAME line that contains no spaces and is short enough to not
               *      exceed the read length.
               *      Conclusion: Just check how BioPython does it:
               *      http://biopython.org/DIST/docs/api/Bio.SeqIO.QualityIO-pysrc.html
               */
              num_reads++;
              line_type = NAME;
            // '@' is a valid quality character
            } else if (line_type == PLUS) {
              line_type = QUAL;
            }
            break;
          case '+':
            if (line_type == SEQ) {
              line_type = PLUS;
            // '+' is a valid quality character
            } else if (line_type == PLUS) {
              line_type = QUAL;
            }
            break;
          default:
            if (line_type == NAME) {
              line_type = SEQ;
            } else if (line_type == PLUS) {
              line_type = QUAL;
            }
        }
        first_char = 0;
      }
    }

    bytes_read = fread(buffer, 1, buffer_size, infile);
  }

  printf("%d\n", num_reads);

  fclose(infile);
  return 0;
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
