#include <iostream>
#include <GASPI.h>
#include <vector>
#include <utility>
#include <ctime>
#include <limits>

#include <skepu>

#include <omp.h>

int main(int argc, char *argv[]){

  int cols;
  int rows;
  size_t iterations;

  if( argc >= 4){
    cols = std::stoi(argv[1]);
    rows = std::stoi(argv[2]);
    iterations = std::stoul(argv[3]);
  }
  else{
    cols = 15;
    rows = 15;
    iterations = 10;
  }


  // 0 args -> we use new with standard args
  // 1 to 2 args -> we use old with standard args
  // 3 args -> we use given args with new
  // 4 or more args -> we use old with given args
  const bool use_old = argc < 4 ?
    argc != 1 :
    argc != 4;

  auto start = std::chrono::system_clock::now();

  skepu::Matrix<float> m1{cols, rows, 1};
  skepu::Matrix<float> m2{cols, rows, 2};
  skepu::Matrix<float> m3{cols, rows, 3};

  const float max = std::numeric_limits<float>::infinity();

  auto mult = skepu::Map<2>([max](float a, float b)  {
    if(a == max || b == max){
        return float{1.01};
      }
    return a * b;
  });

  if(use_old){
    for(int i = 0; i < iterations; i++){
      mult.apply_old(m1, m2, m3);
      mult.apply_old(m2, m1, m3);
      mult.apply_old(m3, m2, m1);
    }
  }
  else{
    for(int i = 0; i < iterations; i++){
      mult.apply_old(m1, m2, m3);
      mult.apply_old(m2, m1, m3);
      mult.apply_old(m3, m2, m1);
    }
  }


  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  float res = m3[0];

  if(use_old){
    std::cout << "Old apply: res = " << res << " calc time = " << elapsed_seconds.count() <<
    " size = " << cols * rows << ", iterations = " << iterations << std::endl;
  }
  else{
    std::cout << "New apply: res = " << res << " calc time = " << elapsed_seconds.count() <<
    " size = " << cols * rows << ", iterations = " << iterations << std::endl;
  }


  return 0;
}
