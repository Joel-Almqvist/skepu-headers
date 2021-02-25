#include <iostream>
#include <skepu>
#include <vector>
#include <utility>


struct Particle
{
	float x, y, z;
	float vx, vy, vz;
	float m;
};



int main(){


  skepu::Matrix<int> m1{4,4,1};
  skepu::Matrix<int> m2{4,4,2};
  skepu::Matrix<int> m3{4,4,3};

	auto mult = skepu::Map<4>([](
		skepu::Index1D i, skepu::Vec<int> a, int b, double c) int {

		/*
		std::cout << "index " << i.i << std::endl;
		std::cout << "proxy " << a.test << std::endl;
		std::cout << "b " << b << std::endl;
		std::cout << "c " << c << std::endl;
		*/

		return 13;
	});

	mult.apply( m1, m2, m3, double{7.77});


  return 0;
}
