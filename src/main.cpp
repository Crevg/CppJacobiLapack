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
int main(int ac, char *av[])
{
  anpi::Matrix<double> A,E;
  A = {{3, 1, 2},{1, 2, 1}, {2, 1, 4}};
 /* A = {
    {4,-30, 60, -35},
    {-30,300, -675, 420},
    {60, -675, 1620, -1050},
    {-35, 420, -1050, 700}
  };*/
  std::vector<double> val;
  anpi::jacobi(A, val, E);
  anpi::eig(A, val, E);

  std::cout << "Eigenvectors: " << std::endl;
  printM(E);
  std::cout << "Eigenvalues: " << std::endl;
  for (size_t i = 0; i < val.size(); ++i){
    std::cout << val[i] << "\t";
  }
  std::cout << "\n";
 
  return 0;
}