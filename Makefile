.POSIX:

.SUFFIXES:

LDLIBS=-lm

.PHONY=all clean

SRC=src/swt.c src/find_seamounts.c
OBJ=$(SRC:.c=.o)

all: swt find_seamounts

swt: src/swt.o
	$(CC) $(LDFLAGS) -o $@ $? $(LDLIBS)

find_seamounts: src/find_seamounts.o
	$(CC) $(LDFLAGS) -o $@ $? $(LDLIBS)

.SUFFIXES: .c .o

.c.o: 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f swt find_seamounts $(OBJ)
