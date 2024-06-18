# define to avoid fully automatic vectorization. If defined, vectorization pragmas will still be vectorized (dependences allowing)
NO_AUTOVEC=-fno-vectorize

CC=clang

#CFLAGS= -march=rv64g -mepi -O2 -g $(NO_AUTOVEC) -mcpu=avispado -mllvm -combiner-store-merging=0 -mllvm -disable-loop-idiom-memset  -mllvm -indvars-widen-indvars=0 -mllvm -disable-loop-idiom-memcpy -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize
CFLAGS= -march=rv64g -mepi -O2 -g $(NO_AUTOVEC) -mcpu=avispado -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize

LDFLAGS = -lm

.PHONY: all 

all: axpy_scalar axpy_intrinsics axpy_autovectorization

axpy_scalar:   main_scalar.o utils.o
	$(CC)  -o $@ $+ $(LDFLAGS)

axpy_intrinsics:   main_vector.o axpy_intrinsic.o utils.o
	$(CC)  -o $@ $+ $(LDFLAGS)

axpy_autovectorization: main_vector.o axpy_autovectorization.o utils.o
	$(CC)  -o $@ $+ $(LDFLAGS)

main_scalar.o: main.c utils.h
	$(CC) $(CFLAGS) -c -o $@ $<

main_vector.o: main.c utils.h
	$(CC) $(CFLAGS) -DVECTOR -c -o $@ $<

axpy_intrinsic.o: axpy.c
	$(CC) $(CFLAGS) -DINTRINSICS -c -o $@ $<

axpy_autovectorization.o: axpy.c
	$(CC) $(CFLAGS) -DAUTOVECTORIZATION -c -o $@ $<

utils.o: utils.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -f *.o
	rm -f axpy_scalar axpy_intrinsics axpy_autovectorization
