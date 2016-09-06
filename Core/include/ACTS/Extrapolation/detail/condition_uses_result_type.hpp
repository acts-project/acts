#ifndef ACTS_CONDITION_USES_RESULT_TYPE_HPP
#define ACTS_CONDITION_USES_RESULT_TYPE_HPP 1

#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  template <typename T, bool has_observer = true>
  struct condition_uses_result_type_impl
  {
    static constexpr bool value = has_result_type_v<observer_type_t<T>>;
  };

  template <typename T>
  struct condition_uses_result_type_impl<T, false> : std::false_type
  {
  };

  template <typename T>
  struct condition_uses_result_type
      : condition_uses_result_type_impl<T, has_observer_type_v<T>>
  {
  };

}  // namespace detail

}  // namespace Acts
#endif  // ACTS_CONDITION_USES_RESULT_TYPE_HPP
