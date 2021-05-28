#ifndef PROXY_MAT_HPP
#define PROXY_MAT_HPP

#include <type_traits>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <deque>
#include <utility>
#include <mutex>

#include <omp.h>
#include <GASPI.h>

#include "matrix.hpp"

namespace skepu{



  template<typename T>
  class Mat{
  private:

    template<int, typename Ret, typename... Func_args>
    friend class Map1D;

    friend class Matrix<T>;

    friend class _gpi::build_tup_util;
    Matrix<T>* owner;


    Mat(Matrix<T>& own){
      owner = &own;
      size = owner->global_size;
    }


  public:
    using is_proxy_type = std::true_type;

    size_t size;

    // Is public to allow for bracket initialization used when creating
    // the arguments tuple in map
    Mat(){}

    T operator()(const size_t row, const size_t col) const {
      return owner->proxy_get(row * owner->col + col);
    }



  };
}
#endif //PROXY_MAT_HPP
