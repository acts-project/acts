#ifndef ACTS_ANY_OF_HPP
#define ACTS_ANY_OF_HPP 1

#include <type_traits>

namespace Acts {

namespace detail {

  // clang-format off
  template <bool... values>
  struct any_of : std::false_type  {};

  template <bool... others>
  struct any_of<true, others...> : public std::true_type {};

  template <bool... others>
  struct any_of<false, others...> : public any_of<others...> {};

  template<bool... values>
  constexpr bool any_of_v = any_of<values...>::value;
  // clang-format on
}  // namespace detail

}  // namespace Acts
#endif  // ACTS_ANY_OF_HPP
