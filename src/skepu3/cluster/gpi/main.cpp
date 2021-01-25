#include <iostream>
#include <GASPI.h>
#include <vector>
#include <utility>

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


  skepu::Matrix<long> m1{4,4,1};
  skepu::Matrix<long> m2{5,5,2};
  skepu::Matrix<long> m3{4,4,3};

  skepu::Matrix<Particle> m4{4,4};


  auto map_add = skepu::Map<2>([](int a, int b) int {
    return a + b;
  });


  return 0;
}
