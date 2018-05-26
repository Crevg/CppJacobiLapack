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


template<typename T>
void eig(const anpi::Matrix<T>& A,
    std::vector<T>& val ,
    anpi::Matrix<T>& E){

    
    /* Locals */
    int N = A.cols();
    int LDA = N;
    int n = N, lda = LDA, info, lwork;
    float wkopt;
    float* work;

    E = A;


    //meter la matriz a en un array para lapack
    float a[LDA*N];
    float w[N];

    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            a[i+j*lda] = A[i][j];
        }
    }

        //guardar en memoria
        lwork = -1;
        ssyev_( "Vectors", "Upper", &n, a, &lda, w, &wkopt, &lwork, &info );
        lwork = (int)wkopt;
        work = (float*)malloc( lwork*sizeof(float) );
        //resolver
        ssyev_( "Vectors", "Upper", &n, a, &lda, w, work, &lwork, &info );
        //revisar convergencia
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        
    //convertir array de lapack a vector de anpi
    for(int i = 0; i < N; ++i){
        val.push_back(w[i]);
    }

    //convertir de nuevo el array de lapack a una matriz de anpi
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            E[i][j] = a[i+j*lda];
        }
    }
        //liberar memoria
        free( (void*)work );
        //exit( 0 );
}


}//anpi

#endif