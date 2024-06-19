CC=clang
CFLAGS=-O3 -ffast-math -I./sdv_tools
VFLAGS=-mepi -mllvm -combiner-store-merging=0 -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize -mcpu=avispado -mllvm -vectorizer-use-vp-strided-load-store -mllvm -enable-mem-access-versioning=0 -mllvm -disable-loop-idiom-memcpy
LIBS=-lm

EXTRAEFLAGS=-DEXTRAE -DEXTRAE_ALWAYS_TRACE -Wl,-rpath -Wl,$(EXTRAE_HOME)/lib -L$(EXTRAE_HOME)/lib -I/$(EXTRAE_HOME)/include -lseqtrace 
VEHAVEFLAGS=-DVEHAVE -I/apps/riscv/vehave/EPI-0.7/development/include/vehave/

all: codes 

test.x: src/test_bl_sgemm.c
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

sgemm.x: src/sgemm.c
	$(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ $(LIBS)

util.x: src/bl_sgemm_util.c
	$(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ $(VEHAVEFLAGS) $(LIBS)

ref.x: src/bl_sgemm_ref.c
	$(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ $(VEHAVEFLAGS) $(LIBS)

%.x: src/%.c
	$(CC) $(CFLAGS) $(VFLAGS) $^ -o $@ $(LIBS)

codes: test.x sgemm.x util.x ref.x

clean:
	rm -f *.x

dist_clean:
	rm -rf *.x extrae_prv_traces vehave_prv_traces *.out *.trace
