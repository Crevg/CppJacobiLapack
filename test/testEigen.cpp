#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <cstdlib>
#include "Matrix.hpp"
#include "Sort.hpp"
#include "RandomSymmetricSqr.hpp"
#include "LapackFunctions.hpp"
#include "Jacobi.hpp"


//create a diagonal matrix from a guiven vector containing diagonal vales
template<typename T>
anpi::Matrix<T> makeDigMatrix(std::vector<T> v) {
	size_t N=v.size();
	anpi::Matrix<T> A(N,N,anpi::DoNotInitialize);//resultant matrix

	for (int i=0;i<int(A.rows());++i){
		for (int j=0;j<int(A.cols());++j){
			if (i==j)
				A[i][j]=v[i];
			else
				A[i][j]=0;
		}
	}
	return A;
}

template<typename T>
anpi::Matrix<T> transpose(anpi::Matrix<T>& E) {
	size_t N=E.rows();
	anpi::Matrix<T> A(N,N,anpi::DoNotInitialize);//transposed matrix
	for(int i =0;i<int(N);i++){
    	for(int j =0;j<int(N);j++){
      	A[j][i]=E[i][j];
    	}
  	}
  	return A;
}





BOOST_AUTO_TEST_SUITE( Eigen )


template<typename T>
void testJacobi() {
	//const T eps = std::numeric_limits<T>::epsilon();
	const T eps=0.5;

	const size_t n=10;
  	anpi::Matrix<T> A,E,D,Et,R;
  	A=anpi::randomSymmetricSqr<T>(n);

  	std::vector<T> v;
  	anpi::jacobi(A,v,E);

  	Et=transpose(E);
  	D=makeDigMatrix(v);
  	R=E*D*Et;
  	//R=transpose(R);
  	for(int i=0;i<int(A.rows());i++){
		for (int j=0;j<int(A.cols());j++){
			std::cout<<"Aij1 "<<A[i][j]<<" Rij "<<R[i][j]<<std::endl;
			BOOST_CHECK( (abs(A[i][j])-abs(R[i][j]))<eps);
		}
	}
}	



template<typename T>
void testLapack() {
	//const T eps = std::numeric_limits<T>::epsilon();
	const T eps = 0.5;
	const size_t n=10;
  	anpi::Matrix<T> A,E,D,Et,R;
  	A=anpi::randomSymmetricSqr<T>(n);

  	std::vector<T> v;
  	anpi::eig(A,v,E);

  	Et=transpose(E);
  	D=makeDigMatrix(v);
  	R=E*D*Et;
  	//R=transpose(R);
  	for(int i=0;i<int(A.rows());i++){
		for (int j=0;j<int(A.cols());j++){
			std::cout<<"Aij2 "<<A[i][j]<<" Rij "<<R[i][j]<<std::endl;
			BOOST_CHECK( (abs(A[i][j])-abs(R[i][j]))<eps );
		}
	}
}	


BOOST_AUTO_TEST_CASE(eigen) {
	
	//testLapack<float>();
	testJacobi<float>();
}

BOOST_AUTO_TEST_SUITE_END()
