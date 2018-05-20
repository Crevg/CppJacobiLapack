#include <cstdlib>
#include "Matrix.hpp"

#ifndef ANPI_LAPACK_HPP
#define ANPI_LAPACK_HPP

namespace anpi
{



template<typename T>
anpi::Matrix<T> randomSymmetricSqr(const size_t N){
    anpi::Matrix<T> M = anpi::Matrix<T>(N, N, 0.0);
    for (int i = 0; i < N; ++i){
        
        for (int k = 0; k < i; ++k){
            M[i][k] = M[k][i];
        }
        for (int j = i; j < N ; ++j){
            M[i][j] = (T) rand() * 100;
        }
    }
    return M;
}

}//anpi

#endif