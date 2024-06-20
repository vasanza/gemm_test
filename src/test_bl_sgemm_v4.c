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
 * test_bl_sgemm.c
 *
 *
 * Purpose:
 * test driver for BLISLAB sgemm routine and reference sgemm routine.
 *
 * Todo:
 *
 *
 * Modification:
 *
 * ----------------------------------------------------------------------
 * gcc test_bl_sgemm_v4.c -o test_bl_sgemm -lm
 *-----------------------------------------------------------------------
 * */


#include "bl_sgemm.h"
#include "bl_config.h"

#define ERROR_TEST

#define TOLERANCE 1E-2
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

void PackWeightLayout(float* dst, const float* src, int nc, int kc, int nr, bool transpose) {
    int index = 0;
    for (int nr_block_start = 0; nr_block_start < nc; nr_block_start += nr) {
        int nr_block_size = nr;
        if((nc - nr_block_start) < nr) nr_block_size = nc - nr_block_start;

        for (int kr_block_start = 0; kr_block_start < kc; kr_block_start++) {
            for (int nr_block_offset = 0; nr_block_offset < nr; nr_block_offset++) {
                if (nr_block_offset >= nr_block_size) {
                    index++;
                    continue;
                }
                int x_idx = transpose ? kr_block_start : (nr_block_start + nr_block_offset);
                int y_idx = transpose ? (nr_block_start + nr_block_offset) : kr_block_start;
                int x_size = transpose ? kc : nc;
                dst[index++] = src[y_idx * x_size + x_idx];
            }
        }
    }
}
//--------------Here is new code from bl_sgemm_util.c bl_sgemm_ref.c--------------
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
/*
 *
 *
 */ 
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

/*
 *
 *
 */
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

// --- Begin Linux build definitions ---

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

// --- End Linux build definitions ---
float bl_clock( void )
{
	return bl_clock_helper();
}
/*------------------------------------------------------------------------------------------
////////////////////////////////////SGEMM CODE from sgemm.c ////////////////////////////////
------------------------------------------------------------------------------------------*/
void AddDot_4x4_opt( int k, float *A, int lda, float *packB, int ldb, float *C, int ldc )
{
   register float C00, C01, C02, C03, C10, C11, C12, C13, C20, C21, C22, C23, C30, C31, C32, C33;
   float *packBp;
   int p;

   C00 = 0.0f;
   C01 = 0.0f;
   C02 = 0.0f;
   C03 = 0.0f;
   C10 = 0.0f;
   C11 = 0.0f;
   C12 = 0.0f;
   C13 = 0.0f;
   C20 = 0.0f;
   C21 = 0.0f;
   C22 = 0.0f;
   C23 = 0.0f;
   C30 = 0.0f;
   C31 = 0.0f;
   C32 = 0.0f;
   C33 = 0.0f;
   for (p = 0; p < k; p++) {
     packBp = &packB[p * 4];

     C00 += A(0, p+0) * packBp[0];
     C01 += A(0, p+0) * packBp[1];
     C02 += A(0, p+0) * packBp[2];
     C03 += A(0, p+0) * packBp[3];
     C10 += A(1, p+0) * packBp[0];
     C11 += A(1, p+0) * packBp[1];
     C12 += A(1, p+0) * packBp[2];
     C13 += A(1, p+0) * packBp[3];
     C20 += A(2, p+0) * packBp[0];
     C21 += A(2, p+0) * packBp[1];
     C22 += A(2, p+0) * packBp[2];
     C23 += A(2, p+0) * packBp[3];
     C30 += A(3, p+0) * packBp[0];
     C31 += A(3, p+0) * packBp[1];
     C32 += A(3, p+0) * packBp[2];
     C33 += A(3, p+0) * packBp[3];

   }
   C(0, 0) += C00;
   C(0, 1) += C01;
   C(0, 2) += C02;
   C(0, 3) += C03;
   C(1, 0) += C10;
   C(1, 1) += C11;
   C(1, 2) += C12;
   C(1, 3) += C13;
   C(2, 0) += C20;
   C(2, 1) += C21;
   C(2, 2) += C22;
   C(2, 3) += C23;
   C(3, 0) += C30;
   C(3, 1) += C31;
   C(3, 2) += C32;
   C(3, 3) += C33;
}


void bl_sgemm_pack(
    int    m,
    int    mr,
    int    n,
    int    nr,
    int    k,
    float *A,
    float *packA,
    int    lda,
    float *B,
    float *packB,
    int    ldb,
    float *C,        // must be aligned
    int    ldc        // ldc must also be aligned
)
{
    int i, j, p;
    int ir, jr;

    // Early return if possible
    if ( m == 0 || n == 0 || k == 0 ) {
        printf( "bl_sgemm(): early return\n" );
        return;
    }
    //------------------------version4------------------------------

    for ( i = 0; i < m; i += DGEMM_MR ) {          // Start 2-nd loop
      for ( j = 0; j < n; j += DGEMM_NR ) {        // Start 1-st loop
           AddDot_4x4_opt( k, &A( i, 0 ), lda, &packB[j * k], ldb, &C( i, j ), ldc );
        }                                          // End   1-st loop
    }                                              // End   2-nd loop
}



/*------------------------------------------------------------------------------------------
//////////////////////////////Main code of matrix multiplication////////////////////////////
------------------------------------------------------------------------------------------*/

void test_bl_sgemm(
        FILE *fp,
        int m,
        int n,
        int k
        ) 
{
    int    i, j, p, nx;
    float *A, *B, *C, *C_ref, *packA, *packB;
    float tmp, error, flops;
    static float ref_beg, ref_time, bl_sgemm_beg, bl_sgemm_time;
    int    nrepeats;
    int    lda, ldb, ldc, ldc_ref;
    float ref_rectime, bl_sgemm_rectime;

    int mr = 4;
    int nr = 4;
    
    A    = (float*)malloc( sizeof(float) * m * k *2);
    B    = (float*)malloc( sizeof(float) * k * n );
    // Allocate packing buffers
    packA  = bl_malloc_aligned( m + mr, k, sizeof(float) );
    packB  = bl_malloc_aligned( k*2, n + nr, sizeof(float) );


    lda = k;
    ldb = n;
    ldc     = n;
    ldc_ref = n;

    C     = bl_malloc_aligned( ldc, n + nr, sizeof(float) );
    // Propósito: La función bl_malloc_aligned asigna un bloque de memoria alineado adecuadamente para operaciones de alta eficiencia.
    // Alineación: Alinear la memoria puede ser crucial para aprovechar ciertas optimizaciones de hardware, especialmente en operaciones con vectores (SIMD) y matrices.
    // Es la dimensión principal o "leading dimension" de la matriz C. En este caso, ldc es igual a m (el número de filas de la matriz C).
    // Esta es la cantidad de columnas a reservar para la matriz C. Se está reservando espacio para n columnas más un adicional de 4 columnas.
    //Este extra podría ser para asegurar un alineamiento correcto en memoria o para evitar problemas de desbordamiento de memoria.
    
    C_ref = (float*)malloc( sizeof(float) * m * n );

    nrepeats = 1; //------------ 3
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
            // A( i, p ) = (float)( i*m + p );
        }
    }
    for ( j = 0; j < n; j ++ ) {
        for ( p = 0; p < k; p ++ ) {
            B( p, j ) = (float)( drand48() );
            // B( p, j ) = (float)( p*n + j );
        }
    }

    for ( j = 0; j < n; j ++ ) {
        for ( i = 0; i < m; i ++ ) {
            C_ref( i, j ) = (float)( 0.0 );	
                C( i, j ) = (float)( 0.0 );	
        }
    }

    PackWeightLayout(packB, B, n, k, nr, false);

    // printf("[B]\n");
    // for(int i = 0; i < k; i++) {
    //   for(int j = 0; j < n; j++) {
    //     printf("%.1f\t", B[i * n + j]);
    //   }
    //   printf("\n");
    // }

    // printf("[packB]\n");
    // for(int i = 0; i < k; i++) {
    //   for(int j = 0; j < nr; j++) {
    //     printf("%.1f\t", packB[i * nr + j]);
    //   }
    //   printf("\n");
    // }

    for ( i = 0; i < nrepeats; i ++ ) {
        bl_sgemm_beg = bl_clock();
        {
            bl_sgemm_pack(
                    m,
                    mr,
                    n,
                    nr,
                    k,
                    A,
                    packA,
                    lda,
                    B,
                    packB,
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
    // printf("ref\n");
    // for(int i = 0; i < m; i++) {
    //   for(int j = 0; j < n; j++) {
    //     printf("%.0f\t", C_ref[i * n + j]);
    //   }
    //   printf("\n");
    // }
    // printf("\n\n");

    // printf("C\n");
    // for(int i = 0; i < m; i++) {
    //   for(int j = 0; j < n; j++) {
    //     printf("%.0f\t", C[i * n + j]);
    //   }
    //   printf("\n");
    // }
    // printf("\n\n");

    // Compute overall floating point operations.
    flops = ( m * n / ( 1000.0 * 1000.0 * 1000.0 ) ) * ( 2 * k );

        // <-------------------------------------------------------------------------printf results real time
    printf( "%5d\t %5d\t %5d\t %5.3lf\t %5.3lf\n", 
            m, n, k, flops / bl_sgemm_rectime, flops / ref_rectime );
    
    // Guardar los resultados en el archivo CSV.
    fprintf(fp, "%d,%d,%d,%5.3lf,%5.3lf\n", m, n, k, flops / bl_sgemm_rectime, flops / ref_rectime);


    free( A     );
    free( packA );
    free( B     );
    free( packB );
    free( C     );
    free( C_ref );
}

/*------------------------------------------------------------------------------------------
/////////////////////////////////////////// main /////////////////////////////////////////
------------------------------------------------------------------------------------------*/
int main( int argc, char *argv[] )
{
    //printf("%%m\t%%n\t%%k\t%%MY_GFLOPS\t%%REF_GFLOPS\n");
    //for(int i = 16; i <= 800; i += 4) {
    //    test_bl_sgemm( i, i, i );
    //}

    //return 0;

    FILE *fp = fopen("results4.csv", "w");
    if (!fp) {
        perror("No se pudo abrir el archivo CSV");
        return 1;
    }

    // Escribir encabezado en el archivo CSV.
    fprintf(fp, "m,n,k,Version4,Version1\n"); //<--------Set actual version

    printf("%%m\t%%n\t%%k\t%%Version4\t%%Version1\n");
    //printf("Start\n");
        for(int i = 16; i < 200; i += 4) {//<---Set max number of iterations and step betwen interations
            test_bl_sgemm(fp,i, i, i);
        }
   

    fclose(fp);
    //printf("ok\n");

    return 0;
}

// gcc test_bl_sgemm_v4.c -o test_bl_sgemm -lm



