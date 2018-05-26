/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <cstdlib>
#include "LapackFunctions.hpp"
#include "RandomSymmetricSqr.hpp"
#include "Jacobi.hpp"
#include "Sort.hpp"
#include <ctime> 

int main(int ac, char *av[])
{
  //declaraciones previas 
  anpi::Matrix<float> A;
  std::vector<float> valJ;
  std::vector<float> valL;
  size_t N = 10;
  std::vector<float> v={1,2,3,4,5,6,7,8,9,10};
  anpi::Matrix<float> EL;
  anpi::Matrix<float> EJ;

  //matriz random
  A = anpi::randomSymmetricSqr<float>(N);

  //ejecución de métodos
  unsigned t0J, t1J;
  t0J=clock();
  anpi::jacobi(A, valJ, EJ);
  t1J=clock();
  double time = (double(t1J-t0J)/CLOCKS_PER_SEC);
  std::cout << "Tiempo Jacobi: " << time << "\n" << std::endl;

  unsigned t0L, t1L;
  t0L=clock();
  anpi::eig(A, valL, EL);
  t1L=clock();
  time = (double(t1L-t0L)/CLOCKS_PER_SEC);
  std::cout << "Tiempo Lapack: " << time << "\n" << std::endl;

  //ordenamiento 
  anpi::sort(valL, EL);
  anpi::sort(valJ, EJ);


  return 0;
}