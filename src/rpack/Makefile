
# Sample settings for GNU compiler + LAPACK
CC=gcc
CFLAGS= -O3
CLIBS=  -llapack -lblas -lm

run_test : test
	./test.x

test: test.c rigid-press.c rigid-press.h
	${CC} ${CFLAGS} *.c -o test.x ${CLIBS}

clean:
	rm -rf *.o *.x


