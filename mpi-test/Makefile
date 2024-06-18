# define to avoid fully automatic vectorization. If defined, vectorization pragmas will still be vectorized (dependences allowing)
NO_AUTOVEC=-fno-vectorize

CC=mpicc

#CFLAGS= -march=rv64g -mepi -O2 -g $(NO_AUTOVEC) -mcpu=avispado -mllvm -combiner-store-merging=0 -mllvm -disable-loop-idiom-memset  -mllvm -indvars-widen-indvars=0 -mllvm -disable-loop-idiom-memcpy -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize
CFLAGS= -march=rv64g -mepi -O2 -g $(NO_AUTOVEC) -mcpu=avispado -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize

LDFLAGS = -lm

.PHONY: all 

all: mpi_test

mpi_test:   main.o utils.o
	OMPI_MPICC=clang $(CC)  -o $@ $+ $(LDFLAGS)

main.o: main.c utils.h
	OMPI_MPICC=clang $(CC) $(CFLAGS) -c -o $@ $<

utils.o: utils.c
	OMPI_MPICC=clang $(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -f *.o
	rm -f mpi_test
