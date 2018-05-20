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
#include "RandomSymmetricSqr.hpp"
#include "Jacobi.hpp"
int main(int ac, char *av[])
{
  anpi::Matrix<double> A,E;
  A = {{1,0,2},{0,2,1},{2,1,1}};
  std::vector<double> val;
  anpi::jacobi(A, val, E);

 
  return 0;
}