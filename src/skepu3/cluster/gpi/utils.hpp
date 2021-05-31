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
}



#endif //UTILS_HPP
