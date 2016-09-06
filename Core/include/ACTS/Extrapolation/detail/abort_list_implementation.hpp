#ifndef ACTS_ABORT_LIST_IMPLEMENTATION_HPP
#define ACTS_ABORT_LIST_IMPLEMENTATION_HPP 1

#include "ACTS/Extrapolation/detail/condition_uses_result_type.hpp"
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  namespace {
    template <bool has_result = true>
    struct condition_caller
    {
      template <typename condition, typename result, typename input>
      static bool
      check(const condition& c, const result& r, input& current)
      {
        typedef observer_type_t<condition>   observer_type;
        typedef result_type_t<observer_type> result_type;

        return c(r.template get<result_type>(), current);
      }
    };

    template <>
    struct condition_caller<false>
    {
      template <typename condition, typename result, typename input>
      static bool
      check(const condition& c, const result&, input& current)
      {
        return c(current);
      }
    };
  }  // end of anonymous namespace

  template <typename... conditions>
  struct abort_list_impl;

  template <typename first, typename... others>
  struct abort_list_impl<first, others...>
  {
    template <typename T, typename input, typename result>
    static bool
    check(const T& conditions_tuple, input& current, const result& r)
    {
      const auto&    this_condition = std::get<first>(conditions_tuple);
      constexpr bool has_result     = condition_uses_result_type<first>::value;
      return condition_caller<has_result>::check(this_condition, r, current)
          || abort_list_impl<others...>::check(conditions_tuple, current, r);
    }
  };

  template <typename last>
  struct abort_list_impl<last>
  {
    template <typename T, typename input, typename result>
    static bool
    check(const T& conditions_tuple, input& current, const result& r)
    {
      constexpr bool has_result     = condition_uses_result_type<last>::value;
      const auto&    this_condition = std::get<last>(conditions_tuple);
      return condition_caller<has_result>::check(this_condition, r, current);
    }
  };

  template <>
  struct abort_list_impl<>
  {
    template <typename T, typename input, typename result>
    static bool
    check(const T&, input&, const result&)
    {
      return false;
    }
  };

}  // namespace

}  // namespace Acts
#endif  // ACTS_ABORT_LIST_IMPLEMENTATION_HPP
