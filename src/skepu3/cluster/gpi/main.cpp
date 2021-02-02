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


  skepu::Matrix<int> m1{4,4,1};
  skepu::Matrix<int> m2{4,4,2};
  skepu::Matrix<int> m3{4,4,3};


  auto add = skepu::Map<2>([](int a, int b) int {
    return a + b;
  });


//	add(m1,m1, m3);
	// TODO compiles with wrong amount of containers
	add.variadic(m1, m2);


	gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
	m1.print();

	//auto map1 = skepu::Map<1413>([](int a, int b) int {
  //  return a + b;
  //});




	//int foo = m1.get(0);
	//printf("Got %d\n", foo);

//gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
  return 0;
}
