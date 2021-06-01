#ifndef PROXY_VEC_HPP
#define PROXY_VEC_HPP

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
  class Vec{
  private:

    template<int, typename Ret, typename... Func_args>
    friend class Map1D;

    friend class Matrix<T>;

    friend class _gpi::build_tup_util;
    Matrix<T>* owner;


    Vec(Matrix<T>& own) :
      size_helper{own.global_size},
      size{size_helper},
      owner{&own} {}

    size_t size_helper;


  public:
    using is_proxy_type = std::true_type;

    const size_t& size;

    // Is public to allow for bracket initialization used when creating
    // the arguments tuple in map
    Vec() :
    size_helper{},
    size{size_helper},
    owner{nullptr}{}


    Vec(const Vec& that) :
      size_helper{that.size_helper},
      size{size_helper},
      owner{that.owner} {}

    Vec<T>& operator=(const Vec<T>&& that){
      owner = that.owner;
      size_helper = that.size_helper;
      return *this;
    }


    T operator[](const size_t i) const {
      return owner->proxy_get(i);
    }

    T operator()(const size_t i) const {
      return owner->proxy_get(i);
    }

  };



  namespace _gpi{

    template<typename T>
    auto is_proxy_t(int) -> decltype(
      std::declval<typename T::is_proxy_type>(),
      std::true_type{}
      ) {
        throw std::logic_error("Function is used only for SFINAE, it should never be called");
      return std::true_type{};
    }

    template<typename T>
    std::false_type is_proxy_t(double){
      throw std::logic_error("Function is used only for SFINAE, it should never be called");
      return std::false_type{};
    }


    template<typename T>
    struct is_skepu_proxy_type{
      using type = decltype(is_proxy_t<T>(int{}));
    };
  }


}
#endif //PROXY_VEC_HPP
