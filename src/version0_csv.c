/*
 * --------------------------------------------------------------------------
 * BLISLAB 
 * --------------------------------------------------------------------------
 * Copyright (C) 2016, The University of Texas at Austin
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name of The University of Texas nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * bl_sgemm.h
 *
 *
 * Purpose:
 * this header file contains all function prototypes.
 *
 * Todo:
 *
 *
 * Modification:
 *
 * 
 * */


//#ifndef BLISLAB_DGEMM_H
//#define BLISLAB_DGEMM_H

// Allow C++ users to include this header file in their source code. However,
// we make the extern "C" conditional on whether we're using a C++ compiler,
// since regular C compilers don't understand the extern "C" construct.
#ifdef __cplusplus
extern "C" {
#endif


//gcc version0.c -o version0 -lm

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
 

#define ERROR_TEST

#define TOLERANCE 1E-2
//-------------------------------------------------------

// Determine the target operating system

#if defined(__linux__)
#define BL_OS_LINUX 1
#else
#error "unsupport OS, this only support Linux"
#endif

// gettimeofday() needs this.
#include <sys/time.h>
#include <time.h>

#define GEMM_SIMD_ALIGN_SIZE 32

#define min( i, j ) ( (i)<(j) ? (i): (j) )

// #define A( i, j )     A[ (j)*lda + (i) ]
// #define B( i, j )     B[ (j)*ldb + (i) ]
// #define C( i, j )     C[ (j)*ldc + (i) ]
// #define C_ref( i, j ) C_ref[ (j)*ldc_ref + (i) ]

#define A( i, j )     A[ (i)*lda + (j) ]
#define B( i, j )     B[ (i)*ldb + (j) ]
#define C( i, j )     C[ (i)*ldc + (j) ]
#define C_ref( i, j ) C_ref[ (i)*ldc_ref + (j) ]

//-------------------------------------------------------------
// from bl_sgemm_ref.c-----
void bl_sgemm_ref(
        int    m,
        int    n,
        int    k,
        float *XA,
        int    lda,
        float *XB,
        int    ldb,
        float *XC,
        int    ldc
        )
{
    // Local variables.
    int    i, j, p;
    float alpha = 1.0, beta = 1.0;

    // Sanity check for early return.
    if ( m == 0 || n == 0 || k == 0 ) return;

    // Reference GEMM implementation.
    for ( i = 0; i < m; i ++ ) {
        for ( p = 0; p < k; p ++ ) {
            for ( j = 0; j < n; j ++ ) {
                XC[ i * ldc + j ] += XA[ i * lda + p ] * XB[ p * ldb + j ];
            }
        }
    }
}

// from bl_sgemm_util.c-----

float *bl_malloc_aligned(
        int    m,
        int    n,
        int    size
        )
{
    float *ptr;
    int    err;

    err = posix_memalign( (void**)&ptr, (size_t)GEMM_SIMD_ALIGN_SIZE, size * m * n );

    if ( err ) {
        printf( "bl_malloc_aligned(): posix_memalign() failures" );
        exit( 1 );    
    }

    return ptr;
}

void bl_sgemm_printmatrix(
        float *A,
        int    lda,
        int    m,
        int    n
        )
{
    int    i, j;
    for ( i = 0; i < m; i ++ ) {
        for ( j = 0; j < n; j ++ ) {
            printf("%lf\t", A[j * lda + i]);
        }
        printf("\n");
    }
}

/*
 * The timer functions are copied directly from BLIS 0.2.0
 *
 */
static float gtod_ref_time_sec = 0.0;

// --- Begin Linux build definitions --

float bl_clock_helper()
{
    float the_time, norm_sec;
    struct timespec ts;

    clock_gettime( CLOCK_MONOTONIC, &ts );

    if ( gtod_ref_time_sec == 0.0 )
        gtod_ref_time_sec = ( float ) ts.tv_sec;

    norm_sec = ( float ) ts.tv_sec - gtod_ref_time_sec;

    the_time = norm_sec + ts.tv_nsec * 1.0e-9;

    return the_time;
}

float bl_clock( void )
{
	return bl_clock_helper();
}

// from sgemm.c----------------------------------
void bl_sgemm(
    int    m,
    int    n,
    int    k,
    float *A,
    int    lda,
    float *B,
    int    ldb,
    float *C,        // must be aligned
    int    ldc        // ldc must also be aligned
)
{
  int    i, j, p;

  // Early return if possible
  if ( m == 0 || n == 0 || k == 0 ) {
    printf( "bl_sgemm(): early return\n" );
    return;
  }
  //---------------------------VERSION 0------------------------------
  for ( i = 0; i < m; i ++ ) {              // Start 2-th loop
      for ( j = 0; j < n; j ++ ) {          // Start 1-nd loop
        for ( p = 0; p < k; p ++ ) {        // Start 0-st loop

              C( i, j ) += A( i, p ) * B( p, j ); //Each operand is a MACRO defined in bl_sgemm() function.

          }                                 // End   0-th loop
      }                                     // End   1-st loop
  }                                         // End   2-nd loop
}
//----------------------------------------------
void computeError(
        int    ldc,
        int    ldc_ref,
        int    m,
        int    n,
        float *C,
        float *C_ref
        )
{
    int    i, j;
    for ( i = 0; i < m; i ++ ) {
        for ( j = 0; j < n; j ++ ) {
            if ( fabs( C( i, j ) - C_ref( i, j ) ) > TOLERANCE ) {
                printf( "C[ %d ][ %d ] != C_ref, %E, %E\n", i, j, C( i, j ), C_ref( i, j ) );
                break;
            }
        }
    }

}

void test_bl_sgemm(
        FILE *fp,
        int m,
        int n,
        int k
        ) 
{
    int    i, j, p, nx;
    float *A, *B, *C, *C_ref;
    float tmp, error, flops;
    float ref_beg, ref_time, bl_sgemm_beg, bl_sgemm_time;
    int    nrepeats;
    int    lda, ldb, ldc, ldc_ref;
    float ref_rectime, bl_sgemm_rectime;

    A    = (float*)malloc( sizeof(float) * m * k );
    B    = (float*)malloc( sizeof(float) * k * n );

    lda = m;
    ldb = k;
    ldc     = m;
    ldc_ref = m;

    C     = bl_malloc_aligned( ldc, n + 4, sizeof(float) );
    // Propósito: La función bl_malloc_aligned asigna un bloque de memoria alineado adecuadamente para operaciones de alta eficiencia.
    // Alineación: Alinear la memoria puede ser crucial para aprovechar ciertas optimizaciones de hardware, especialmente en operaciones con vectores (SIMD) y matrices.
    // Es la dimensión principal o "leading dimension" de la matriz C. En este caso, ldc es igual a m (el número de filas de la matriz C).
    // Esta es la cantidad de columnas a reservar para la matriz C. Se está reservando espacio para n columnas más un adicional de 4 columnas.
    //Este extra podría ser para asegurar un alineamiento correcto en memoria o para evitar problemas de desbordamiento de memoria.
    
    C_ref = (float*)malloc( sizeof(float) * m * n );

    nrepeats = 5; //------------ 3
    // Medición del Rendimiento: Repetir las operaciones varias veces ayuda a obtener una medida más precisa del rendimiento.
    // Las mediciones pueden variar debido a factores como la carga del sistema, caché hits/misses, etc. 
    // Al repetir la operación y tomar el tiempo más bajo (bl_sgemm_rectime y ref_rectime), se intenta obtener una estimación más fiable del rendimiento óptimo.
    // Promedio o Mejor Tiempo: El código actual toma el menor tiempo registrado en las repeticiones para cada operación (bl_sgemm y bl_sgemm_ref).
    // Esto se hace para minimizar la influencia de cualquier anomalía temporal que pueda ocurrir durante una de las ejecuciones (como un breve aumento en la carga del sistema).

    srand48 (time(NULL));

    // Randonly generate points in [ 0, 1 ].
    for ( p = 0; p < k; p ++ ) {
        for ( i = 0; i < m; i ++ ) {
            A( i, p ) = (float)( drand48() );	
        }
    }
    for ( j = 0; j < n; j ++ ) {
        for ( p = 0; p < k; p ++ ) {
            B( p, j ) = (float)( drand48() );
        }
    }

    for ( j = 0; j < n; j ++ ) {
        for ( i = 0; i < m; i ++ ) {
            C_ref( i, j ) = (float)( 0.0 );	
                C( i, j ) = (float)( 0.0 );	
        }
    }

    for ( i = 0; i < nrepeats; i ++ ) {
        bl_sgemm_beg = bl_clock();
        {
            bl_sgemm(
                    m,
                    n,
                    k,
                    A,
                    lda,
                    B,
                    ldb,
                    C,
                    ldc
                    );
        }
        bl_sgemm_time = bl_clock() - bl_sgemm_beg;

        if ( i == 0 ) {
            bl_sgemm_rectime = bl_sgemm_time;
        } else {
            bl_sgemm_rectime = bl_sgemm_time < bl_sgemm_rectime ? bl_sgemm_time : bl_sgemm_rectime;
        }
    }

#ifdef ERROR_TEST
    for ( i = 0; i < nrepeats; i ++ ) {
        ref_beg = bl_clock();
        {
            bl_sgemm_ref(
                    m,
                    n,
                    k,
                    A,
                    lda,
                    B,
                    ldb,
                    C_ref,
                    ldc_ref
                    );
        }
        ref_time = bl_clock() - ref_beg;

        if ( i == 0 ) {
            ref_rectime = ref_time;
        } else {
            ref_rectime = ref_time < ref_rectime ? ref_time : ref_rectime;
        }
    }

    computeError(
            ldc,
            ldc_ref,
            m,
            n,
            C,
            C_ref
            );
#endif

    // Compute overall floating point operations.
    flops = ( m * n / ( 1000.0 * 1000.0 * 1000.0 ) ) * ( 2 * k );
    
    // <-------------------------------------------------------------------------printf results real time
    printf( "%5d\t %5d\t %5d\t %5.3lf\t %5.3lf\n", 
            m, n, k, flops / bl_sgemm_rectime, flops / ref_rectime );
    
    // Guardar los resultados en el archivo CSV.
    fprintf(fp, "%d,%d,%d,%5.3lf,%5.3lf\n", m, n, k, flops / bl_sgemm_rectime, flops / ref_rectime);


    free( A     );
    free( B     );
    free( C     );
    free( C_ref );
}

int main( int argc, char *argv[] )
{
    //printf("%%m\t%%n\t%%k\t%%MY_GFLOPS\t%%REF_GFLOPS\n");
    //for(int i = 16; i <= 800; i += 4) {
    //    test_bl_sgemm( i, i, i );
    //}

    //return 0;

    FILE *fp = fopen("results.csv", "w");
    if (!fp) {
        perror("No se pudo abrir el archivo CSV");
        return 1;
    }

    // Escribir encabezado en el archivo CSV.
    fprintf(fp, "m,n,k,Version0,Version1\n"); //<----------------------------Set actual version

    printf("%%m\t%%n\t%%k\t%%Version0\t%%Version1\n");
    //printf("Start\n");
    //<--------------------------------------------------------------------------------
    //for(int i = 16; i <= 800; i += 4) {
    
    // Prueba 1: A(m=4k,K=4k) * B(K=4k,n=4k) = C(m=4k,n=4k)
    //for (int j = 10; j <= 800; j+= 100){
        //for(int i = 16; i <= 800; i += 4) {
        for(int i = 20; i < 100; i += 2) {//<----------------------------------Set max number of iterations and step betwen interations
            test_bl_sgemm(fp,i, i, i);
            //test_bl_sgemm(i, i, i);
        }
    //}

    // Prueba 2: A(m=4k,K=4k) * B(K=4k,n=11k) = C(m=4k,n=11k)
    //for(int i = 16; i <= 4000; i += 4) {
    //    test_bl_sgemm(fp,i, i, 11000);
    //}

    fclose(fp);
    //printf("ok\n");

    return 0;
}


