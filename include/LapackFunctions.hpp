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
extern void ssyev( char* jobz, char* uplo, int* n, float* a, int* lda,
                float* w, float* work, int* lwork, int* info );


extern void print_matrix( char* desc, int m, int n, float* a, int lda );


template<typename T>
void eig(const anpi::Matrix<T>& A,
    std::vector<T>& val ,
    anpi::Matrix<T>& E){

    /*int N, LDA = A.cols();
    int LWORK = 3*N-1;
    int INFO;
    anpi::Matrix<T> U = getUpper(A);
    std::vector<T> W = std::vector<T>(N);
    std::vector<T> WORK = std::vector<T>(LWORK);

    ssyev ('V', 'U', N, U, LDA, W, WORK, LWORK, INFO);*/

    /* Locals */
    int N = 5;
    int LDA = N;
        int n = N, lda = LDA, info, lwork;
        float wkopt;
        float* work;
        /* Local arrays */
        float w[N];
        float a[LDA*N] = {
            1.96f,  0.00f,  0.00f,  0.00f,  0.00f,
           -6.49f,  3.80f,  0.00f,  0.00f,  0.00f,
           -0.47f, -6.39f,  4.17f,  0.00f,  0.00f,
           -7.20f,  1.50f, -1.51f,  5.70f,  0.00f,
           -0.65f, -6.34f,  2.67f,  1.80f, -7.10f
        };
        /* Executable statements */
        printf( " SSYEV Example Program Results\n" );
        /* Query and allocate the optimal workspace */
        lwork = -1;
        ssyev( "Vectors", "Upper", &n, a, &lda, w, &wkopt, &lwork, &info );
        lwork = (int)wkopt;
        work = (float*)malloc( lwork*sizeof(float) );
        /* Solve eigenproblem */
        ssyev( "Vectors", "Upper", &n, a, &lda, w, work, &lwork, &info );
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

//sacar la matriz superior
template<typename T>
anpi::Matrix<T> getUpper(const anpi::Matrix<T>& A){
    int t = A.cols();
    anpi::Matrix<T> U = anpi::Matrix<T>(t, t, 0.0);
    
    for(int j = 0; j < t; ++j){
         
        for(int i = j; i < t; ++i){
            U[i][j] = A[i][j];
        }
    }
    return U;
}


}//anpi

#endif