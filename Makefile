CC = gcc
CFLAGS = -Wall

all:
	$(CC) $(CFLAGS) readsfq.c -o readsfq
	$(CC) $(CFLAGS) -shared -fPIC swalign.c -o libswalign.so -lm
.PHONY: all

clean:
	rm -f readsfq libswalign.so
.PHONY: clean

