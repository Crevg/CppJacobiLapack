#include <cstdlib>
#include <math.h>
#include "Matrix.hpp"
#include <stdlib.h>
#include <stdio.h>

#ifndef ANPI_LAPACK_HPP
#define ANPI_LAPACK_HPP


namespace anpi
{

//Calcula los eigenvalores y vectores, la escogimos porque sirve para calcular los eigenvalores y vectores
//de una matriz real simétrica y porque los parámetros son los mas convenientes para nuestro propósito
// como la definición de la matriz superior o inferior, que en otro método se definia como un vector;
// pero en este se define como una matriz, y por lo tanto es más fácil manejar 
//std::ssyev ( const char*, const char*, const int* , const std::vector<T>*, const int*, const std::vector<T>*,  const std::vector<T>*, const int*, const int*);

/* SSYEV prototype */

extern "C"
{
	extern void ssyev_( char* jobz, char* uplo, int* n, float* a, int* lda,
                float* w, float* work, int* lwork, int* info );
}

extern void print_matrix( char* desc, int m, int n, float* a, int lda );


template<typename T>
void eig(const anpi::Matrix<T>& A,
    std::vector<T>& val ,
    anpi::Matrix<T>& E){

    /* Locals */
    int N = 10;
    int LDA = N;
    int n = N, lda = LDA, info, lwork;
    float wkopt;
    float* work;

    //meter la matriz a en un array para lapack
    float a[LDA*N];
    int cont = 0;

    int col = A.cols();
    for(int i = 0; i < col; ++i){
        for(int j = 0; j < col; ++j){
            a[cont] = A[i][j];
            ++cont;
        }
    }
        /* Local arrays */
        float w[N];

        /* Executable statements */
        printf( " SSYEV Example Program Results\n" );
        /* Query and allocate the optimal workspace */
        lwork = -1;
        ssyev_( "Vectors", "Upper", &n, a, &lda, w, &wkopt, &lwork, &info );
        lwork = (int)wkopt;
        work = (float*)malloc( lwork*sizeof(float) );
        /* Solve eigenproblem */
        ssyev_( "Vectors", "Upper", &n, a, &lda, w, work, &lwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        /* Print eigenvalues */
        print_matrix( "Eigenvalues", 1, n, w, 1 );
        /* Print eigenvectors */
        print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
        /* Free workspace */
        free( (void*)work );
        exit( 0 );
} /* End of SSYEV Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, float* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
                printf( "\n" );
        }

    
}


}//anpi

#endif