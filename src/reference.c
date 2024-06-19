#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
 
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

void test_bl_sgemm(
        //FILE *fp,
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
    //fprintf(fp, "%d,%d,%d,%5.3lf,%5.3lf\n", m, n, k, flops / bl_sgemm_rectime, flops / ref_rectime);


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

    //FILE *fp = fopen("results.csv", "w");
    //if (!fp) {
    //    perror("No se pudo abrir el archivo CSV");
    //    return 1;
    //}

    // Escribir encabezado en el archivo CSV.
    //fprintf(fp, "m,n,k,Version0,Version1\n"); //<----------------------------Set actual version

    printf("%%m\t%%n\t%%k\t%%Version0\t%%Version1\n");
    //printf("Start\n");
    //<--------------------------------------------------------------------------------
    //for(int i = 16; i <= 800; i += 4) {
    
    // Prueba 1: A(m=4k,K=4k) * B(K=4k,n=4k) = C(m=4k,n=4k)
    //for (int j = 10; j <= 800; j+= 100){
        //for(int i = 16; i <= 800; i += 4) {
        for(int i = 20; i < 1000; i += 20) {//<----------------------------------Set max number of iterations and step betwen interations
            //test_bl_sgemm(fp,i, i, i);
            test_bl_sgemm(i, i, i);
        }
    //}

    // Prueba 2: A(m=4k,K=4k) * B(K=4k,n=11k) = C(m=4k,n=11k)
    //for(int i = 16; i <= 4000; i += 4) {
    //    test_bl_sgemm(fp,i, i, 11000);
    //}

    //fclose(fp);
    //printf("ok\n");

    return 0;
}


