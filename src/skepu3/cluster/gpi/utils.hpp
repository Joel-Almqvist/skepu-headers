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
  };

}


namespace skepu::_gpi{

  // Returns a tuple type containing N copies of T
  // Credits for these structs goes to rahnema1 from Stackoverflow

  // Usage:
  // tuple_of<3, int> tup = make_tuple(1,1,1);
  template <size_t N,typename T>
  struct tuple_n{
      template< typename...Args>
  		using type = typename tuple_n<N-1, T>::template type<T, Args...>;
  };

  template <typename T>
  struct tuple_n<0, T> {
      template<typename...Args>
  		using type = std::tuple<Args...>;
  };

  template <size_t N,typename T>
  using tuple_of = typename tuple_n<N, T>::template type<>;



  // The dummy is used to allow compile time evaluation of the
  // bool not_done
  template <int ctr, bool done>
  struct dummy;

  template <int ctr>
  struct dummy<ctr, true> {

      template <typename Func, typename...Args, typename... Exp>
      static auto exec(Func func, std::tuple<Args...>& tup, Exp&... exp)
      -> typename std::remove_reference<decltype(std::get<0>(tup))>::type
       {

        const bool not_done = ctr < sizeof...(Args) - 1;
        dummy<ctr + 1, not_done>::exec(func, tup, exp..., std::get<ctr>(tup));
      }
  };

  template <int ctr>
  struct dummy<ctr, false> {

      template<typename Func, typename Tup, typename...Exp>
      static auto exec(Func func, Tup& tup, Exp&... exp)
      -> typename std::remove_reference<decltype(std::get<0>(tup))>::type
      {
        return func(exp...);
        }
  };

/* Should be struct based for compile time usage
  template<typename First, typename... Rest>
  bool is_skepu(First& first, Rest&... rest){
    return First::is_skepu_container && is_skepu(rest...);
  }

  template<typename Last>
  bool is_skepu(Last& last){
    return Last::is_skepu_container;
  }
*/


// **********************
template<typename T>
auto get_val_t(int) -> decltype(
  std::declval<typename T::is_skepu_container>(),
  (typename T::value_type){}
  ) {
  return (typename T::value_type){};
}

template<typename T>
T get_val_t(double){
  return T{};
}

template<typename T>
struct get_skepu_val_t{
  //using type = decltype(std::declval<get_val_t<T>(0)>());
  using type = decltype(get_val_t<T>(int{}));
};
// ***************************

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
  //using type = decltype(std::declval<get_val_t<T>(0)>());
  using type = decltype(is_skepu_helper<T>(int{}));
};

// ******************************



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


// Given a lambda Func with n + 1 arguments create a tuple type
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
