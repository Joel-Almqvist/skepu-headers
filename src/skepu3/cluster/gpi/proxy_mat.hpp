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


    Mat(Matrix<T>& own) :
      cols_helper{own.col},
      rows_helper{own.row},
      size_helper{own.global_size},
      cols{cols_helper},
      rows{rows_helper},
      size{size_helper},
      owner{&own} {}


    size_t cols_helper;
    size_t rows_helper;
    size_t size_helper;

  public:
    using is_proxy_type = std::true_type;

    // The user may only access non readable values
    const size_t& size;
    const size_t& cols;
    const size_t& rows;

    // Is public to allow for bracket initialization used when creating
    // the arguments tuple in map
    Mat() :
    cols_helper{},
    rows_helper{},
    size_helper{},
    cols{cols_helper},
    rows{rows_helper},
    size{size_helper},
    owner{nullptr} {}


    Mat(const Mat& that) :
      cols_helper{that.cols_helper},
      rows_helper{that.rows_helper},
      size_helper{that.size_helper},
      cols{cols_helper},
      rows{rows_helper},
      size{size_helper},
      owner{that.owner} {}

    Mat<T>& operator=(const Mat<T>&& that){
      owner = that.owner;
      size_helper = that.size_helper;
      cols_helper = that.cols_helper;
      rows_helper = that.rows_helper;

      return *this;
    }

    // 2D accessing pattern
    T operator()(const size_t row, const size_t col) const {
      return owner->proxy_get(row * owner->col + col);
    }

    // 1D accessing pattern
    T operator[](const size_t i) const {
      return owner->proxy_get(i);
    }



  };
}
#endif //PROXY_MAT_HPP
