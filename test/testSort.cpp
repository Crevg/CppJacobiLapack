#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <cstdlib>
#include "Matrix.hpp"
#include "Sort.hpp"




template<typename T>
void testSort() {

	const T eps = std::numeric_limits<T>::epsilon();

	

	anpi::Matrix<T> original,sorted;
	std::vector<T> values;
	values={8,4,6,2,1,9,7};
	original = {{22,33,44,55,66,77,88},
				{23,34,45,56,67,78,89}, 
				{24,35,46,57,68,79,81},
				{25,36,47,58,69,71,82},
				{26,37,48,59,61,72,83},
				{27,38,49,51,62,73,84},
				{28,39,41,52,63,74,85}};

	sorted =   {{77,22,88,44,33,55,66},
				{78,23,89,45,34,56,67}, 
				{79,24,81,46,35,57,68},
				{71,25,82,47,36,58,69},
				{72,26,83,48,37,59,61},
				{73,27,84,49,38,51,62},
				{74,28,85,41,39,52,63}};

	anpi::sort(values, original);

	for(int i=0;i<int(original.rows());i++){
		for (int j=0;j<int(original.cols());j++){
			BOOST_CHECK( abs(original[i][j]-sorted[i][j])<eps );
		}
	}




}	

BOOST_AUTO_TEST_SUITE( Sort )

BOOST_AUTO_TEST_CASE(sort) {
	testSort<float>();
}

BOOST_AUTO_TEST_SUITE_END()
