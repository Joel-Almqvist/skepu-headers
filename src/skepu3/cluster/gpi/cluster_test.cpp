#include <iostream>
#include <GASPI.h>
#include <vector>
#include <utility>

#include <matrix.hpp>
#include <reduce.hpp>
#include <map.hpp>
#include <omp.h>


int main(){

  auto square = skepu::Map<1>([](long a) long {
    return a * a;
  });


  auto add = skepu::Map<2>([](long a, long b) long {
    return a + b;
  });


  auto mult = skepu::Map<2>([](long a, long b) long {
    return a * b;
  });

  // Constructor: Matrix(M, N initial value)
  // The Matrix interface will need to be tweaked to match SkePU's
  skepu::Matrix<long> m1{1000, 100, 1};
  skepu::Matrix<long> m2{1000, 100, 2};
  skepu::Matrix<long> m3{1000, 100, 3};

  for(int i = 0; i < 10; i ++){
  add(m1, m2, m3);
  add(m2, m1, m3);
  add(m3, m2, m1);
}

  int res = m3.get(0);
  std::cout << "Res = " << res << std::endl;


  return 0;
}
