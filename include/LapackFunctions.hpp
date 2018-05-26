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


/**
   *Declaración de la función ssyev e Lapack para 
   *calcular los eigenvalores y eigenvectores 
   *
   * @param[in] JOBZ = 'N':  calcula solo eigenvalores;
   *                 = 'V':  calcula eigenvalores y eigenvectores. 
   * @param[in] IPLO = 'U':  Triangular superior;
   *                 = 'L':  tringular inferior. 
   * @param[in] N orden de la matriz 
   * @param[in,out] A  arreglo real de tamaño (LDA, N)
   * @param[in] LDA dimensiones del arreglo A.  LDA >= max(1,N). 
   * @param[out] W un arreglo real, tamaño (N)
   * @param[out] WORK arreglro real, tamaño (MAX(1,LWORK))
   * @param[out] LWORK es un int LWORK >= max(1,3*N-1).
   *             si LWORK = -1, Calcula el valor mas óptimo para WORK
   * @param[out]INFO es un int
                = 0:  successful exit   
                < 0:  if INFO = -i, el i-ésimo argumento tiene un valor no permitido
                > 0:  if INFO = i, el algoritmo no converge
   */
extern "C"
{
	extern void ssyev_( char* jobz, char* uplo, int* n, float* a, int* lda,
                float* w, float* work, int* lwork, int* info );
}

/**
   *Calculo los eigenvalores y eigenvectores por medio del método de lapack
   *
   * @param[in] A matriz cuadrada real 
   * @param[out] vector con los eigenvalores de A
   * @param[out] Matriz con los eigenvectores de A
   */
template<typename T>
void eig(const anpi::Matrix<T>& A,
    std::vector<T>& val ,
    anpi::Matrix<T>& E){

    
    /* Declaraciones */
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