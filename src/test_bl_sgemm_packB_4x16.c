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
 * gcc test_bl_sgemm_packB_4x16.c sgemm.c bl_sgemm_util.c bl_sgemm_ref.c RvvSgemm4x16.S -o test_bl_sgemm -lm
 * ./test_bl_sgemm
 * python3 plot_results.py
 *-----------------------------------------------------------------------
 * 
 * */



#include "bl_sgemm.h"

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

void test_bl_sgemm(
        //FILE *fp,
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
    int nr = 16;
    
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

    nrepeats = 1;//------------------------- 3
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
    //fprintf(fp, "%d,%d,%d,%5.3lf,%5.3lf\n", m, n, k, flops / bl_sgemm_rectime, flops / ref_rectime);

    free( A     );
    free( packA );
    free( B     );
    free( packB );
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
    fprintf(fp, "m,n,k,Version5,Version1\n"); //<----------------------------Set actual version

    //printf("%%m\t%%n\t%%k\t%%Version0\t%%Version1\n");
    printf("Start\n");
    //<--------------------------------------------------------------------------------
    //for(int i = 16; i <= 800; i += 4) {
    
    // Prueba 1: A(m=4k,K=4k) * B(K=4k,n=4k) = C(m=4k,n=4k)
    //for (int j = 10; j <= 800; j+= 100){
        //for(int i = 16; i <= 800; i += 4) {
        for(int i = 20; i < 100; i += 20) {//<----------------------------------Set max number of iterations and step betwen interations
            test_bl_sgemm(fp,i, i, i);
        }
    //}

    // Prueba 2: A(m=4k,K=4k) * B(K=4k,n=11k) = C(m=4k,n=11k)
    //for(int i = 16; i <= 4000; i += 4) {
    //    test_bl_sgemm(fp,i, i, 11000);
    //}

    fclose(fp);
    printf("ok\n");

    return 0;
}

