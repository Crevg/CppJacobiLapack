#include <cstdlib>
#include <algorithm> 
#include "Matrix.hpp"
#include <time.h>
#ifndef ANPI_SORT_HPP
#define ANPI_SORT_HPP

namespace anpi
{

template<typename T>
void swapCols(int a,int b,anpi::Matrix<T>& E){
    int i;
    for(i=0;i<int(E.rows());i++){
        T temp=E[i][b];
        E[i][b]=E[i][a];
        E[i][a]=temp;
    }
}

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