#include <cstdlib>
#include <algorithm> 
#include "Matrix.hpp"
#include <time.h>
#ifndef ANPI_SORT_HPP
#define ANPI_SORT_HPP

namespace anpi
{

/**
  *Intercambia las columnas indicadas en los parametros a y b de la 
  *matriz E elemento por elemento.
  *@param[in] a primera columna por ser intercambiada
  *@param[in] b columna para ser intercambiada con a
  *@param[out] E matriz resultante de intercambiar las columnas
  */

template<typename T>
void swapCols(int a,int b,anpi::Matrix<T>& E){
    int i;
    for(i=0;i<int(E.rows());i++){
        T temp=E[i][b];
        E[i][b]=E[i][a];
        E[i][a]=temp;
    }
}



/**
  *Reordenar el vector de eigenvalores en orden decreciente. 
  *Al mismo tiempo realiza intercambios en las columnas de 
  *la matriz de eigenvectores segun los cambios en el orden
  *del vector. El ordenamiento del vector se realiza mediante 
  *bubble sort. Las columnas de la matriz 
  *@param[out] val vector de eigenvalores
  *@param[out] E matriz con eigenvectores
  */
template<typename T>
void sort(std::vector<T>&  val, anpi::Matrix<T>& E){
   int i, j;
   for (i = 0; i < int(val.size())-1; i++)      
 
       // Last i elements are already in place   
       for (j = 0; j < int(val.size())-i-1; j++) 
           if (abs(val[j]) < abs(val[j+1])){
                T temp =val[j];
                val[j]=val[j+1];
                val[j+1]=temp;
                swapCols(j,j+1,E);
           }
}//sort
 


}//anpi

#endif