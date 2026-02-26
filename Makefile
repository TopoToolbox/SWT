.POSIX:

.SUFFIXES: .c .o .out .chk

LDLIBS=-lm

.PHONY=all clean check

SRC=src/swt.c src/find_seamounts.c
OBJ=$(SRC:.c=.o)

all: swt find_seamounts

swt: src/swt.o
	$(CC) $(LDFLAGS) -o $@ $? $(LDLIBS)

find_seamounts: src/find_seamounts.o
	$(CC) $(LDFLAGS) -o $@ $? $(LDLIBS)

# This is the standard C compilation rule
.c.o: 
	$(CC) $(CFLAGS) -c $< -o $@

# Testing
TEST=test/block.txt test/coeffcomb.txt test/coeffint.txt test/coeffraw.txt test/count.txt test/wavfilt.txt

TESTOUT=$(TEST:.txt=.out)
TESTCHK=$(TEST:.txt=.chk)

check: $(TESTCHK)

$(TESTOUT): swt test/input.txt
	./swt test/input.txt test/coeffraw.out test/coeffint.out test/wavfilt.out test/coeffcomb.out test/block.out test/count.out > test/output.txt

# .chk is a pseudo-suffix: nothing ever gets generated with this
# suffix. Instead we use it to trigger this rule that compares the
# stored output ($*.txt) with the generated output ($<).
.out.chk:
	diff $*.txt $<

clean:
	rm -f swt find_seamounts $(OBJ) $(TESTOUT) look.temp lookB.temp
