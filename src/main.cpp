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
 


  //pruebas para reconstruir A como E*D*Et con jacobi
 size_t NE=A.rows();
  anpi::Matrix<float> EJt(NE,NE,anpi::DoNotInitialize);
  for(int i =0;i<int(NE);i++){
      for(int j =0;j<int(NE);j++){
        EJt[j][i]=EJ[i][j];
      }
    }

  anpi::Matrix<float> Diag(NE,NE,anpi::DoNotInitialize);//resultant matrix

  for (int i=0;i<int(EJ.rows());++i){
    for (int j=0;j<int(EJ.cols());++j){
      if (i==j)
        Diag[i][j]=valJ[i];
      else
        Diag[i][j]=0;
    }
  }
  
  anpi::Matrix<float> R;
  R=EJ*Diag*EJt;
  std::cout << "A ES: " << "\n" << std::endl;
  printM(A);
  std::cout << "\n";
  std::cout << "EJ ES: " << "\n" << std::endl;
  printM(EJ);
  std::cout << "\n";

  std::cout << "EJt ES: " << "\n" << std::endl;
  printM(EJt);
  std::cout << "\n";
  
  std::cout << "Diag ES: " << "\n" << std::endl;
  printM(Diag);
  std::cout << "\n";
  std::cout << "Resultado de EJ*Diag*EJt ES: " << "\n" << std::endl;
  printM(R);
  std::cout << "\n";


   //pruebas para reconstruir A como E*D*Et con lapack
  anpi::Matrix<float> ELt(NE,NE,anpi::DoNotInitialize);
  for(int i =0;i<int(NE);i++){
      for(int j =0;j<int(NE);j++){
        ELt[j][i]=EL[i][j];
      }
    }

  anpi::Matrix<float> Diag2(NE,NE,anpi::DoNotInitialize);//resultant matrix

  for (int i=0;i<int(EL.rows());++i){
    for (int j=0;j<int(EL.cols());++j){
      if (i==j)
        Diag2[i][j]=valL[i];
      else
        Diag2[i][j]=0;
    }
  }
  
  anpi::Matrix<float> R2;
  R2=EL*Diag2*ELt;
  std::cout << "A ES: " << "\n" << std::endl;
  printM(A);
  std::cout << "\n";
  std::cout << "EL ES: " << "\n" << std::endl;
  printM(EL);
  std::cout << "\n";

  std::cout << "ELt ES: " << "\n" << std::endl;
  printM(ELt);
  std::cout << "\n";
  
  std::cout << "Diag ES: " << "\n" << std::endl;
  printM(Diag);
  std::cout << "\n";
  std::cout << "Resultado de EL*Diag2*ELt ES: " << "\n" << std::endl;
  printM(R2);
  std::cout << "\n";

  return 0;
}