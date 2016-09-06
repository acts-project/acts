#ifndef ACTS_ABORT_LIST_CREATOR_HPP
#define ACTS_ABORT_LIST_CREATOR_HPP 1

#include <tuple>
#include <type_traits>
#include "ACTS/Extrapolation/AbortConditions.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Extrapolation/detail/Extendable.hpp"
#include "ACTS/Extrapolation/detail/condition_uses_result_type.hpp"
#include "ACTS/Utilities/detail/MPL/boost_mpl_helper.hpp"
#include "ACTS/Utilities/detail/MPL/has_duplicates.hpp"

namespace Acts {

namespace {
  template <typename T,
            typename input,
            typename result,
            typename = decltype(std::declval<T>().
                                operator()(std::declval<const result&>(),
                                           std::declval<input&>()))>
  std::true_type
  test_condition_with_result(int);

  template <typename, typename, typename>
  std::false_type test_condition_with_result(...);

  template <typename T,
            typename input,
            typename = decltype(std::declval<T>().
                                operator()(std::declval<input&>()))>
  std::true_type
  test_condition_without_result(int);

  template <typename, typename>
  std::false_type test_condition_without_result(...);

  // clang-format off
  template <typename T, typename input, bool has_result = false>
  struct conditions_traits_check_impl : decltype(test_condition_without_result<T, input>(0)) {};

  template <typename T, typename input>
  struct conditions_traits_check_impl<T, input, true>
      : decltype(test_condition_with_result<T,input,
          typename detail::result_type_t<detail::observer_type_t<T>>>(0))
  {
  };

  template <typename T, typename input>
  struct conditions_traits_checker
      : conditions_traits_check_impl<T, input, detail::condition_uses_result_type<T>::value>
  {
  };
  // clang-format on

  template <bool has_result = true>
  struct ConditionCaller
  {
    template <typename condition, typename result, typename input>
    static bool
    check(const condition& c, const result& r, input& current)
    {
      typedef detail::observer_type_t<condition>   observer_type;
      typedef detail::result_type_t<observer_type> result_type;

      return c(r.template get<result_type>(), current);
    }
  };

  template <>
  struct ConditionCaller<false>
  {
    template <typename condition, typename result, typename input>
    static bool
    check(const condition& c, const result&, input& current)
    {
      return c(current);
    }
  };

  template <typename... conditions>
  struct AbortListImpl;

  template <typename first, typename... others>
  struct AbortListImpl<first, others...>
  {
    template <typename T, typename input, typename result>
    static bool
    check(const T& conditions_tuple, input& current, const result& r)
    {
      const auto&    this_condition = std::get<first>(conditions_tuple);
      constexpr bool has_result
          = detail::condition_uses_result_type<first>::value;
      return ConditionCaller<has_result>::check(this_condition, r, current)
          || AbortListImpl<others...>::check(conditions_tuple, current, r);
    }
  };

  template <typename last>
  struct AbortListImpl<last>
  {
    template <typename T, typename input, typename result>
    static bool
    check(const T& conditions_tuple, input& current, const result& r)
    {
      constexpr bool has_result
          = detail::condition_uses_result_type<last>::value;
      const auto& this_condition = std::get<last>(conditions_tuple);
      return ConditionCaller<has_result>::check(this_condition, r, current);
    }
  };

  template <>
  struct AbortListImpl<>
  {
    template <typename T, typename input, typename result>
    static bool
    check(const T&, input&, const result&)
    {
      return false;
    }
  };
}

template <typename... conditions>
struct AbortList : private detail::Extendable<conditions...>
{
private:
  static_assert(not detail::has_duplicates_v<conditions...>,
                "same abort conditions type specified several times");

  typedef detail::type_collector_t<detail::observer_type_extractor,
                                   conditions...>
      observers;

  using detail::Extendable<conditions...>::tuple;

public:
  typedef detail::unpack_boost_set_as_tparams_t<ObserverList, observers>
      observer_list_type;
  using detail::Extendable<conditions...>::get;

  template <typename input, typename result_t>
  bool
  operator()(input& current, const result_t& r) const
  {
    static_assert(
        all_of_v<conditions_traits_checker<conditions, input>::value...>,
        "not all abort conditions support the specified input");

    AbortListImpl<conditions...>::check(tuple(), current, r);
  }
};

}  // namespace Acts

#endif  // ACTS_ABORT_LIST_CREATOR_HPP
