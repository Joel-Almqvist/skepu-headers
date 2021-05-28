#ifndef UTILS_HPP
#define UTILS_HPP

#include <type_traits>
#include <tuple>


namespace skepu{

  // Does this type trait exists elsewhere? Do we need to have it here?
  template<typename T>
  struct is_skepu_container : std::false_type {};

  struct Index1D{
    size_t i;

    using is_skepu_index = std::true_type;
    using is_skepu_1D_index = std::true_type;
  };

  struct Index2D{
    size_t row;
    size_t col;

    using is_skepu_index = std::true_type;
    using is_skepu_2D_index = std::true_type;
  };

}


namespace skepu::_gpi{

template<typename T>
auto is_skepu_helper(int) -> decltype(
  std::declval<typename T::is_skepu_container>())
  {
    return std::true_type{};
  }

template<typename T>
std::false_type is_skepu_helper(double){
  return std::false_type{};
}

template<typename T>
struct is_skepu{
  using type = decltype(is_skepu_helper<T>(int{}));
};



// This type trait was written by Potatoswatter from Stackoverflow:
template< typename t, std::size_t n, typename = void >
struct function_argument_type;

template< typename r, typename ... a, std::size_t n >
struct function_argument_type< r (*)( a ... ), n >
    { typedef typename std::tuple_element< n, std::tuple< a ... > >::type type; };

template< typename r, typename c, typename ... a, std::size_t n >
struct function_argument_type< r (c::*)( a ... ), n >
    : function_argument_type< r (*)( a ... ), n > {};

template< typename r, typename c, typename ... a, std::size_t n >
struct function_argument_type< r (c::*)( a ... ) const, n >
    : function_argument_type< r (c::*)( a ... ), n > {};

template< typename ftor, std::size_t n >
struct function_argument_type< ftor, n,
    typename std::conditional< false, decltype( & ftor::operator () ), void >::type >
    : function_argument_type< decltype( & ftor::operator () ), n > {};

// endof Potatoswatter's code


// Given a lambda or functor with n + 1 arguments create a tuple type
// containing the arguments.
template<typename Func, int n, typename ... Rest>
struct get_tup_t{

  using type = typename get_tup_t<Func, n - 1,
    typename function_argument_type<Func, n >::type, Rest...>::type;
};


  template<typename Func, typename ... Rest>
  struct get_tup_t<Func, 0, Rest...>{

    using type = std::tuple<typename function_argument_type<Func, 0 >::type,
      Rest...>;
  };
}



#endif //UTILS_HPP
