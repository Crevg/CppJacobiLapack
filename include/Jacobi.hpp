#include <cstdlib>
#include "Matrix.hpp"
#include <cmath>

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

template<typename T>
std::vector<size_t> getMax(const anpi::Matrix<T>&A){
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
                if (abs(A(i,j)) > abs(A(maxI, maxJ))){
                    maxI = i;
                    maxJ = j;
                }
            }
        }
    }
    std::vector<size_t> v;
    v.push_back(maxI);
    v.push_back(maxJ);
    return v; 
}

template<typename T>
void jacobi(const anpi::Matrix<T>& A, std::vector<T> &val, anpi::Matrix<T> &E){
    T t ,ttemp , tau, c, s;
    std::vector<size_t> max = anpi::getMax(A);
    size_t row = max[0];
    size_t col = max[1];
    std::cout << "Posicion a cambiar: (" << row << ", " << col <<
        ") y (" << col << ", " << row << ")\n";
    if (A[row][col] == 0){
        std::cout << "c'est finit" << std::endl;
        return;
    }
    tau = (A[col][col] - A[row][row])/(2*A[row][col]);
    t = -tau - sqrt(1+ tau*tau);
    ttemp = -tau + sqrt(1+ tau*tau);
    if (abs(ttemp) < abs(t)){
        t = ttemp;
    }
    c = 1/ sqrt(1+t*t);
    s = c * t;
    std::cout << "Metadatos: " << std::endl;
    std::cout << "T: " << t << "\ttau: " << tau <<
        "\tc: " << c << "\ts: " << s << "\t(Ttemp: " << ttemp << ")\n"; 
    size_t N = A.rows();
    anpi::Matrix<T> Pi = anpi::Matrix<T>(N,N,0.0);
    anpi::Matrix<T> PiT = anpi::Matrix<T>(N,N,0.0);
    for (int i = 0; i < N; ++i){
        Pi(i,i) = PiT(i,i) = 1;
    }
    Pi(row, row) = Pi(col, col) = c;
    PiT(row, row) = PiT(col, col) = c;
    Pi(row, col) = PiT(col, row) = s;
    Pi(col, row) = PiT(row, col) = -s;

    std::cout << "Matrices P: " << std::endl;
    std::cout << "Pi: " << std::endl;
    printM(Pi);
    std::cout << "PiT: " << std::endl;
    printM(PiT);

    E = PiT*A*Pi;
    std::cout << "E: " << std::endl;
    printM(E);
}

}//anpi

#endif