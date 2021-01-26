#include <iostream>
#include <GASPI.h>
#include <vector>
#include <utility>

#include <matrix.hpp>
#include <reduce.hpp>
#include <map.hpp>
#include <omp.h>
#include <ctime>

int main(int argc, char *argv[]){

  

  auto start = std::chrono::system_clock::now();
  

  auto square = skepu::Map<1>([](long a)  {
    return a * a;
  });


  auto add = skepu::Map<2>([](long a, long b)  {
    return a + b;
  });


  auto mult = skepu::Map<2>([](long a, long b)  {
    return a * b;
  });

  int i;
   i = atoi(argv[1]);
  //i = 32;

  // Constructor: Matrix(M, N initial value)
  // The Matrix interface will need to be tweaked to match SkePU's
  skepu::Matrix<long> m1{10000, 1000, 1};
  skepu::Matrix<long> m2{10000, 1000, 1};
  skepu::Matrix<long> m3{10000, 1000, 1};

  for(; i > 0; i --){
  add(m1, m2, m3);
  add(m2, m1, m3);
  add(m3, m2, m1);
}

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  //int res = m3.get(0);
  int res = 0;
  std::cout << "Res = " << res << " calc time = " << elapsed_seconds.count() << std::endl;


  return 0;
}
