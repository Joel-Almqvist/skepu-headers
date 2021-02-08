#include <iostream>
#include <GASPI.h>
#include <vector>
#include <utility>

#include<environment.hpp>
#include <matrix.hpp>
#include <reduce.hpp>
#include <map.hpp>
#include <omp.h>
//#include <filter.hpp>

struct Particle
{
	float x, y, z;
	float vx, vy, vz;
	float m;
};


int main(){

  skepu::Matrix<int> m1{4,4,1};
  skepu::Matrix<int> m2{4,4,2};
  skepu::Matrix<int> m3{5,5,3};
	skepu::Matrix<int> m4{4,4,4};
	skepu::Matrix<int> m5{6,3,5};

	auto square = skepu::Map<1>([](int a) int {
		return a * a;
	});

  auto add2 = skepu::Map<2>([](int a, int b) int {
    return a + b;
  });

	auto add3 = skepu::Map<3>([](int a, int b, int c) int {
    return a + b + c;
  });


	auto add4 = skepu::Map<4>([](int a, int b, int c, int d) int {
		return a + b + c + d;
	});


	add2(m1, m2, m2);
	add3(m1, m1, m2, m3);
	add3(m1, m1, m2, m3);
	add3(m1, m1, m2, m3);
	add3(m1, m1, m2, m3);
	add3(m1, m1, m2, m3);
	add3(m1, m1, m2, m3);

	int foo = m2.get(4);
	foo = m2.get(5);
	foo = m2.get(10);


	m1.print();

  return 0;
}
