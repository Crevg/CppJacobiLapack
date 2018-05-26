#include <cstdlib>
#include "Matrix.hpp"
#include <cmath>
#include <limits>

#ifndef ANPI_JACOBI_HPP
#define ANPI_JACOBI_HPP

namespace anpi
{

/***
 * Verifica que el par ordenado (row, col) no existan respectivamente en
 * los vectores iVec y jVec y retorna un booleano indicandolo.
 * @param [in] row numero a verificar en el vector de filas
 * @param [in] col numero a verificar en el vector de columnas
 * @param [in] iVec vector de filas
 * @param [in] jVec vector de columnas
 * 
 * @return false si existen como par, true si no
 * */
bool notInV(size_t row, size_t col, const std::vector<size_t> &iVec, const std::vector<size_t> &jVec){
    size_t n = iVec.size();
    for (size_t i = 0; i<n; ++i){
        if (iVec[i] == row && jVec[i] == col){
            return false;
        }else if (jVec[i] == row && iVec[i] == col){
            return false;
        }
    }
    return true;
}


/***
 * Obtiene el par (fila, columna) que contengan el dato de mayor
 * valor absoluto afuera de la diagonal
 * @param [in] A matriz con los valores a verificar
 * @param [out] rows vector donde se agrega la fila del par ordenado encontrado
 * @param [out] cols vector donde se agrega la columna del par ordenado encontrado
 * 
 * @throws anpi::Exception si la matriz no es cuadrada.
 * */
template<typename T>
void getMax(const anpi::Matrix<T>&A, std::vector<size_t>& rows, 
                            std::vector<size_t> &cols){
    size_t maxI = 0;
    size_t maxJ= 1;
    size_t N = A.rows();
    size_t M = A.cols();
    
    if (N != M){
        throw anpi::Exception("Matrix has to be nxn");
    }
    for (size_t i = 0; i < N; ++i){
        for (size_t j = 0; j < N; ++j){
            if (i != j){
                if (abs(A(i,j)) > abs(A(maxI, maxJ)) && notInV(i, j, rows, cols)){
                    maxI = i;
                    maxJ = j;
                }
            }
        }
    }
    rows.push_back(maxI);
    cols.push_back(maxJ);
}


/***
 * Obtiene la matriz de eigenvectores y los eigenvalores de una matriz 
 * por medio del metodo de transformaciones de Jacobi
 * @param [in] A matriz con los valores a verificar
 * @param [out] val vector con los eigenvalores
 * @param [out] E matriz de eigenvectores (en las columnas)
 * */
template<typename T>
void jacobi(const anpi::Matrix<T>& A, std::vector<T> &val, anpi::Matrix<T> &E){
    /* declara las variables necesarias */
    T t ,ttemp , tau, c, s, S, eps;
    std::vector<size_t> rows, cols;
    size_t row, col, N, vLen, sweeps;
    anpi::Matrix<T> Ap;
    S = 1;
    T Sprev = 0;
    N = A.rows();
    Ap = A;
    /* inicializa E con la matriz identidad NxN*/
    E = anpi::Matrix<T>(N,N,0.0);
    for (size_t i = 0 ; i < N; ++i){
        E(i,i) = 1.0;
    }

    eps = sqrt(std::numeric_limits<T>::epsilon());
    sweeps = 1;
    int n = 0;
    while(S > eps){    
          
        /*calcula posiciones a cambiar */
        anpi::getMax(Ap, rows, cols);
        vLen = rows.size();
        row = rows[vLen-1]; //last element
        col = cols[vLen-1]; //last element
        tau = (Ap[col][col] - Ap[row][row])/(2*Ap[row][col]);
        t = -tau - sqrt(1+ tau*tau);
        ttemp = -tau + sqrt(1+ tau*tau);
        if (abs(ttemp) < abs(t)){
            t = ttemp;
        }
        c = 1/ sqrt(1+t*t);
        s = c * t;

        /* crear y llena la matriz Pi y su transpuesta */ 

        anpi::Matrix<T> Pi = anpi::Matrix<T>(N,N,0.0);
        anpi::Matrix<T> PiT = anpi::Matrix<T>(N,N,0.0);

        for (size_t j = 0; j < N; ++j){
            Pi(j,j) = PiT(j,j) = 1;
        }
        Pi(row, row) = Pi(col, col) = c;
        PiT(row, row) = PiT(col, col) = c;
        Pi(row, col) = PiT(col, row) = s;
        Pi(col, row) = PiT(row, col) = -s;

        /* crear la matriz de la siguiente iteración */ 
        Ap = PiT*Ap*Pi;
        
        /*criterio de convergencia */ 
        S = 0;
        for (size_t l = 0; l < N; ++l){
            for (size_t k = 0; k < N; ++k){
                if (l != k){
                    S += abs(Ap[l][k]*Ap[l][k]);
                }
            }
        }
        S -= 2*abs(Ap[row][col]*Ap[row][col]);
        S = abs(S);
        
        E = E*Pi; //guarda los eigenvectores en las columnas

        /*check if sweeps is finished */
        if (vLen == N*(N-1)/2){
            sweeps+=1;
            rows.clear();
            cols.clear();

            /* si 50 barridos seguidos mantiene el mismo 
             * criterio de convergencia entonces termine
             */
            if (S == Sprev){
                n++;
                if (n == 50){
                    break;
                }
            }else{
                n = 0;
            }
            Sprev = S;
        }    
    }
    for (size_t i = 0; i < N; ++i){
        val.push_back(Ap(i,i));
    }
}

}//anpi

#endif