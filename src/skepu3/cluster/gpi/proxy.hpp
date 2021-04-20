#ifndef PROXY_HPP
#define PROXY_HPP

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


    Vec(Matrix<T>& own){
      owner = &own;
      size = owner->global_size;
    }


  public:
    using is_proxy_type = std::true_type;

    size_t size;

    // Is public to allow for bracket initialization used when creating
    // the arguments tuple in map
    Vec(){
    }

    T operator[](const size_t i){
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





    // Functions to easier access the proxy type. Not currently used.
    template<typename T>
    auto make_proxy_t(int) -> decltype(
      std::declval<typename T::is_skepu_container>(),
      Vec<T>{}
      ) {
        throw std::logic_error("Function is used only for SFINAE, it should never be called");
      return Vec<T>{};
    }

    template<typename T>
    T make_proxy_t(double){
      throw std::logic_error("Function is used only for SFINAE, it should never be called");
      return T{};
    }

    // Takes a Skepu container type and returns a proxy type
    template<typename T>
    struct get_skepu_proxy_t{
      using type = decltype(make_proxy_t<T>(int{}));
    };




    // Takes in a list of types and returns a tuple type of said elements with
    // every skepu container replaced with a proxy type.
    /* Ex:
    using foo_t = typename _gpi::pack<std::tuple<>, DestCont, Conts...>::type;

    bool const b = std::is_same<foo_t, std::tuple<_gpi::placeholder,
    _gpi::placeholder, _gpi::placeholder>>::value;
    */
    template<typename ...Args>
    struct pack{};

    template<typename ... Unpacked, typename Curr, typename ...Packed>
    struct pack<std::tuple<Unpacked...>, Curr, Packed...>{
      using helper = typename std::tuple<Unpacked...,
        typename get_skepu_proxy_t<Curr>::type>;
      using type = typename pack<helper, Packed...>::type;
    };

    template<typename ... Unpacked, typename Curr>
    struct pack<std::tuple<Unpacked...>, Curr>{
      using type = std::tuple<Unpacked..., typename get_skepu_proxy_t<Curr>::type>;
    };
  }


}
#endif //PROXY_HPP
