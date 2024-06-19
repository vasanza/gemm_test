CC=clang
CFLAGS=-O3 -ffast-math -I./sdv_tools
VFLAGS=-mepi -mllvm -combiner-store-merging=0 -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize -mcpu=avispado -mllvm -vectorizer-use-vp-strided-load-store -mllvm -enable-mem-access-versioning=0 -mllvm -disable-loop-idiom-memcpy
LIBS=-lm

EXTRAEFLAGS=-DEXTRAE -DEXTRAE_ALWAYS_TRACE -Wl,-rpath -Wl,$(EXTRAE_HOME)/lib -L$(EXTRAE_HOME)/lib -I/$(EXTRAE_HOME)/include -lseqtrace 
VEHAVEFLAGS=-DVEHAVE -I/apps/riscv/vehave/EPI-0.7/development/include/vehave/

all: codes 

reference.x: src/my_sgemm.c
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

reference-i-extrae.x: src/my_sgemm.c
	$(CC) $(CFLAGS) $^ -o $@ $(EXTRAEFLAGS) $(LIBS)

reference-vec.x: src/my_sgemm.c
	$(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ $(LIBS)

reference-vec-i-vehave.x: src/my_sgemm.c
	$(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ $(VEHAVEFLAGS) $(LIBS)

%-i-vehave.x: src/%.c
	$(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ $(VEHAVEFLAGS) $(LIBS)

%-i-extrae.x: src/%.c
	$(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ $(EXTRAEFLAGS) $(LIBS)

%.x: src/%.c
	$(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ $(LIBS)

codes: reference.x reference-vec.x increase-vec.x increase-vl.x flex-datatype.x
codes-vehave: reference-vec-i-vehave.x light-i-vehave.x increase-vec-i-vehave.x increase-vl-i-vehave.x flex-datatype-i-vehave.x
codes-extrae: reference-i-extrae.x increase-vec-i-extrae.x increase-vl-i-extrae.x flex-datatype-i-extrae.x friendly-i-extrae.x

clean:
	rm -f *.x

dist_clean:
	rm -rf *.x extrae_prv_traces vehave_prv_traces *.out *.trace
