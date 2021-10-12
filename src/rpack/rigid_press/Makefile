
# Sample settings for GNU compiler + LAPACK
CC=gcc
CFLAGS= -O3 -Wall -DROPT_DEBUG
CLIBS=  -llapack -lblas -lm

run_test : test
	./test.x

test: test.c rigid_press.c rigid_press.h symmetrization.c spglib
	${CC} ${CFLAGS} test.c rigid_press.c symmetrization.c d_algebra.c  *.o -o test.x ${CLIBS}

spglib:
	${CC} ${CFLAGS} -c ../spglib_src/*.c  ${CLIBS}

clean:
	rm -rf *.o *.x __pycache__ build  _rigid*


