#ifndef UTILS_HPP
#define UTILS_HPP

#include <type_traits>


namespace skepu{

  // Does this type trait exists elsewhere? Do we need to have it here?
  template<typename T>
  struct is_skepu_container : std::false_type {};

}

#endif //UTILS_HPP
