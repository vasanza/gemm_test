#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>

#include "utils.h" 


void axpy_intrinsics(double a, double *dx, double *dy, int n);

// Ref version
void axpy_ref(double a, double *dx, double *dy, int n) {
   int i;
   for (i=0; i<n; i++) {
      dy[i] += a*dx[i];
   }
}

void init(double *pv, long n, double value)
{
   for (int i=0; i<n; i++) pv[i]= value;
}

int main(int argc, char *argv[])
{
    double a=1.0;
    long n;

    if (argc == 2)
       n = atol(argv[1]); // input argument: vector size in Ks
    else{
       printf ("USAGE: axpy num_elements\n");
       return 1;
     }

    /* Allocate the source and result vectors */
    double *dx     = (double*)malloc(n*sizeof(double));
    double *dy     = (double*)malloc(n*sizeof(double));
    double *dy_ref = (double*)malloc(n*sizeof(double));

    uint64_t start_cycles, start_instret, end_cycles, end_instret;

    init(dx, n, 1.0);
    init(dy, n, 2.0);
    printf ("doing reference axpy\n");
    __asm__ __volatile__("rdinstret %0" : "=r"(start_instret));
    __asm__ __volatile__("rdcycle %0" : "=r"(start_cycles));
    axpy_ref(a, dx, dy, n); 
    __asm__ __volatile__("rdcycle %0" : "=r"(end_cycles));
    __asm__ __volatile__("rdinstret %0" : "=r"(end_instret));
    uint64_t cycles = end_cycles - start_cycles;
    uint64_t instret = end_instret - start_instret;
    printf ("reference: done\n");
    printf ("reference: cycles=%lu instructions=%lu\n", cycles, instret);
    capture_ref_result(dy, dy_ref, n);
#ifdef VECTOR
    init(dx, n, 1.0);
    init(dy, n, 2.0);
    printf ("doing vector axpy\n");
    __asm__ __volatile__("rdinstret %0" : "=r"(start_instret));
    __asm__ __volatile__("rdcycle %0" : "=r"(start_cycles));
    axpy_vector(a, dx, dy, n);
    __asm__ __volatile__("rdcycle %0" : "=r"(end_cycles));
    __asm__ __volatile__("rdinstret %0" : "=r"(end_instret));
    cycles = end_cycles - start_cycles;
    instret = end_instret - start_instret;
    printf ("vector: done\n");
    printf ("vector: cycles=%lu instructions=%lu\n", cycles, instret);
    test_result(dy, dy_ref, n);
#endif
    free(dx); free(dy); free(dy_ref);
    return 0;
}
