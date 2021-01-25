#include <iostream>
#include <GASPI.h>
#include <vector>
#include <utility>

#include <matrix.hpp>
#include <reduce.hpp>
#include <map.hpp>
#include <omp.h>




int main(){


  auto square = skepu::Map<1>([](int a) int {
    return a * a;
  });


  auto add = skepu::Map<2>([](int a, int b) int {
    return a + b;
  });


  auto mult = skepu::Map<2>([](int a, int b) int {
    return a * b;
  });



	int run_test = 0;

	if(run_test == 0){
		// Test repeated use of the same container

		// Constructor: Matrix(M, N initial value)
		skepu::Matrix<int> m2{4,4,2};
		skepu::Matrix<int> m22{4,4,2};

		// The Matrix interface will need to be tweaked to match SkePU's

		// 2 * 2
		square(m2, m2);
		// 4 * 4
		square(m2, m2);
		// 16 * 16
		square(m2, m2);
		// 256 * 256 = 65536
		square(m2, m2);

		mult(m22, m22, m22);
		mult(m22, m22, m22);
		mult(m22, m22, m22);
		mult(m22, m22, m22);

		m2.print();
		m22.print();
	}


	else if(run_test == 1){
		// Different individual values test

		skepu::Matrix<int> m1{4,4,1};
		skepu::Matrix<int> m3{5,5,3};
		skepu::Matrix<int> m4{4,4,4};


		for(int i = 0; i < 16; i++){
			m1.set(i, i);
		}

		add(m1, m1, m3);
		m1.print();


	}

	else{
		// Combined operations with different sizes of matrices

		// Note that m3 is larger than m1 and m4 and hence can never be the
		// destination of a map


		skepu::Matrix<int> m1{4,4,1};
		skepu::Matrix<int> m3{5,5,3};
		skepu::Matrix<int> m4{4,4,4};

		// m1 = m4 + m4
		add(m1, m4, m4);

		// m4 = m1^2
		square(m4, m1);

		// m4 = m4 * m3
		mult(m4, m4, m3);

		// m1 = m4 + m3
		add(m1, m4, m3);

		// 195
		m1.print();
	}


  return 0;
}
