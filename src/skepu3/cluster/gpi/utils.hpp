#ifndef UTILS_HPP
#define UTILS_HPP

#include <type_traits>
#include<tuple>

namespace skepu{

  // Does this type trait exists elsewhere? Do we need to have it here?
  template<typename T>
  struct is_skepu_container : std::false_type {};
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
      static auto exec(Func func, std::tuple<Args...> tup, Exp... exp)
      -> decltype(std::get<0>(tup))
       {

        const bool not_done = ctr < sizeof...(Args) - 1;
        dummy<ctr + 1, not_done>::exec(func, tup, exp..., std::get<ctr>(tup));
      }
  };

  template <int ctr>
  struct dummy<ctr, false> {

      template<typename Func, typename Tup, typename...Exp>
      static auto exec(Func func, Tup tup, Exp... exp)
      -> decltype(std::get<0>(tup))
      {
        func(exp...);
        }
  };

}



#endif //UTILS_HPP
