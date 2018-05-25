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
  //A = {{3, 1, 2},{1, 2, 1}, {2, 1, 4}};
 /*A = {
    {4,-30, 60, -35},
    {-30,300, -675, 420},
    {60, -675, 1620, -1050},
    {-35, 420, -1050, 700}
  };*/
  std::vector<float> valJ = {0,0,0,0,0,0,0,0,0,0};
  std::vector<float> valL = {0,0,0,0,0,0,0,0,0,0};
  size_t N = 10;
  std::vector<float> v={1,2,3,4,5,6,7,8,9,10};
  anpi::Matrix<float> EL(N,N,anpi::DoNotInitialize);
  anpi::Matrix<float> EJ(N,N,anpi::DoNotInitialize);

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

  //impresión jacobi
  std::cout << "Eigenvectors Jacobi: " << "\n" << std::endl;
  printM(EJ);
  std::cout << "\n";
  std::cout << "Eigenvalues Jacobi: " << "\n" << std::endl;
  for (size_t i = 0; i < valJ.size(); ++i){
    std::cout << valJ[i] << "   ";
  }
  std::cout << "\n";
  std::cout << "\n";

  //impresión Lapack
  std::cout << "Eigenvectors Lapack: " << "\n" << std::endl;
  printM(EL);
  std::cout << "\n";
  std::cout << "Eigenvalues Lapack: " << "\n" << std::endl;
  for (size_t i = 0; i < valL.size(); ++i){
    std::cout << valL[i] << "   ";
  }
  std::cout << "\n";
  std::cout << "\n";
 
  return 0;
}