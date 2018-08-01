CFLAGS = -Wall

all:
	gcc $(CFLAGS) readsfq.c -o readsfq
.PHONY: all

clean:
	rm -f readsfq
.PHONY: clean

