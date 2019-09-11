CFLAGS = -Wall

all:
	gcc $(CFLAGS) readsfq.c -o readsfq
	gcc $(CFLAGS) -shared -fPIC swalign.c -o libswalign.so -lm
.PHONY: all

clean:
	rm -f readsfq libswalign.so
.PHONY: clean

