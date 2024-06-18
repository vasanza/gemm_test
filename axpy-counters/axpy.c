#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#if __riscv_vector_version==800
#define BROADCAST_f64 __builtin_epi_vfmv_v_f_1xf64
#else
#define BROADCAST_f64 __builtin_epi_vbroadcast_1xf64
#endif

#ifdef AUTOVECTORIZATION
void axpy_vector(double a, double *dx, double *dy, int n) {
   int i;
   #pragma clang loop vectorize(enable)
   for (i=0; i<n; i++) {
      dy[i] += a*dx[i];
   }
}

#endif

#ifdef INTRINSICS
void axpy_vector(double a, double *dx, double *dy, int n) {
  int i;


  long gvl = __builtin_epi_vsetvl(n, __epi_e64, __epi_m1);
  __epi_1xf64 v_a = BROADCAST_f64(a, gvl);
  
  for (i = 0; i < n;) {
    gvl = __builtin_epi_vsetvl(n - i, __epi_e64, __epi_m1);
    __epi_1xf64 v_dx = __builtin_epi_vload_1xf64(dx, gvl);
    dx+=gvl;
    __epi_1xf64 v_dy = __builtin_epi_vload_1xf64(dy, gvl);
    __epi_1xf64 v_res = __builtin_epi_vfmacc_1xf64(v_dy, v_a, v_dx, gvl);
    __builtin_epi_vstore_1xf64(dy, v_res, gvl);
    dy+=gvl;
    i += gvl;
  }
}
#endif
