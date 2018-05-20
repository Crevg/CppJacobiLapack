#include <cstdlib>
#include "Matrix.hpp"
#include <time.h>
#ifndef ANPI_RANDOM_SYMMETRIC_SQR_HPP
#define ANPI_RANDOM_SYMMETRIC_SQR_HPP

namespace anpi
{



template<typename T>
anpi::Matrix<T> randomSymmetricSqr(const size_t N){
    srand(time(NULL));
    T low = -1000;
    T high = 1000;
    anpi::Matrix<T> M = anpi::Matrix<T>(N, N, 0.0);
    for (size_t i = 0; i < N; ++i){
        
        for (size_t k = 0; k < i; ++k){
            M[i][k] = M[k][i];
        }
        for (size_t j = i; j < N ; ++j){
            M[i][j] = low + static_cast <T> (rand()) /(static_cast <T> (RAND_MAX/(high-low)));
        }
    }
    return M;
}

}//anpi

#endif