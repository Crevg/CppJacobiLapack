#include <cstdlib>
#include "Matrix.hpp"
#include <cmath>
#include <limits>

#ifndef ANPI_JACOBI_HPP
#define ANPI_JACOBI_HPP

namespace anpi
{

template <typename T>
void printM(const anpi::Matrix<T> &A){
    size_t n = A.rows();
    for (size_t i = 0; i < n; ++i){
        for (size_t j = 0; j < n; ++j){
            std::cout << A(i,j) << "\t";
        }
        std::cout << std::endl;
    }
}
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
                /*std::cout << "i,j: " << abs(A(i,j)) <<
                " iMax, jMax curr: "  << abs(A(maxI, maxJ)) <<
                " comparison " << (abs(A(i,j)) > abs(A(maxI, maxJ))) <<
                " notInV? " << notInV(maxI, maxJ, rows, cols) << std::endl;
                */
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

template<typename T>
void jacobi(const anpi::Matrix<T>& A, std::vector<T> &val, anpi::Matrix<T> &E){
    T t ,ttemp , tau, c, s;
    std::vector<size_t> rows, cols;
    size_t row, col, N, vLen, sweeps;
    T S= 1;
    N = A.rows();
    E = A;
    size_t iterator = N*(N-1)/2;
    T eps = sqrt(std::numeric_limits<T>::epsilon());
    sweeps = 1;
    while(S > eps){
        
        /*calcula posiciones a cambiar */
        anpi::getMax(E, rows, cols);
        vLen = rows.size();
        row = rows[vLen-1]; //last element
        col = cols[vLen-1]; //last element
        std::cout << "Posicion a cambiar: (" << row << ", " << col <<
            ") y (" << col << ", " << row << ")\n";

        /*calcula los parametros */
        tau = (E[col][col] - E[row][row])/(2*E[row][col]);
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
        E = PiT*E*Pi;
        
        /*criterio de convergencia */ 
        S = 0;
        for (size_t l = 0; l < N; ++l){
            for (size_t k = 0; k < N; ++k){
                if (l != k){
                    S += E[l][k]*E[l][k];
                }
            }
        }
        S -= 2*E[row][col]*E[row][col];
        S = abs(S);
        std::cout << "S' igual: " << S << std::endl; 
        
        /*check if sweeps is finished */
        if (vLen == N*(N-1)/2){
            sweeps+=1;
            rows.clear();
            cols.clear();
        }
    }

    printM(E);


}

}//anpi

#endif