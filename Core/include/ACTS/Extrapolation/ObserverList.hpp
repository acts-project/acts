#ifndef ACTS_OBSERVER_LIST_HPP
#define ACTS_OBSERVER_LIST_HPP 1

#include <iostream>
#include <tuple>
#include <utility>
#include "ACTS/Extrapolation/detail/Extendable.hpp"
#include "ACTS/Utilities/detail/MPL/all_of.hpp"
#include "ACTS/Utilities/detail/MPL/boost_mpl_helper.hpp"
#include "ACTS/Utilities/detail/MPL/has_duplicates.hpp"
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace {
  using detail::all_of_v;
  using detail::has_duplicates_v;
  using detail::type_collector_t;
  using detail::result_type_extractor;
  using detail::unpack_boost_set_as_tparams_t;

  template <typename T,
            typename input,
            typename result,
            typename = decltype(std::declval<T>().
                                operator()(std::declval<const input&>(),
                                           std::declval<const input&>(),
                                           std::declval<result&>()))>
  std::true_type
  test_with_result(int);

  template <typename, typename, typename>
  std::false_type test_with_result(...);

  template <typename T,
            typename input,
            typename = decltype(std::declval<T>().
                                operator()(std::declval<const input&>(),
                                           std::declval<const input&>()))>
  std::true_type
  test_without_result(int);

  template <typename, typename>
  std::false_type test_without_result(...);

  // clang-format off
  template <typename T, typename input, bool has_result = false>
  struct observer_traits_check_impl : decltype(test_without_result<T, input>(0)) {};
  // clang-format on

  template <typename T, typename input>
  struct observer_traits_check_impl<T, input, true>
      : decltype(
            test_with_result<T, input, typename result_type_extractor::type<T>>(
                0))
  {
  };

  template <typename T, typename input>
  struct observer_traits_checker
      : observer_traits_check_impl<T,
                                   input,
                                   result_type_extractor::has_type<T>::value>
  {
  };

  template <bool has_result = true>
  struct ObserverCaller
  {
    template <typename observer, typename result, typename input>
    static void
    observe(const observer& obs,
            const input&    current,
            const input&    previous,
            result&         r)
    {
      obs(current,
          previous,
          r.template get<typename result_type_extractor::type<observer>>());
    }
  };

  template <>
  struct ObserverCaller<false>
  {
    template <typename observer, typename result, typename input>
    static void
    observe(const observer& obs,
            const input&    current,
            const input&    previous,
            result&)
    {
      obs(current, previous);
    }
  };

  template <typename... observers>
  struct ObserverListImpl;

  template <typename first, typename... others>
  struct ObserverListImpl<first, others...>
  {
    template <typename T, typename result, typename input>
    static void
    observe(const T&     obs_tuple,
            const input& current,
            const input& previous,
            result&      r)
    {
      const auto& this_observer = std::get<first>(obs_tuple);
      ObserverCaller<result_type_extractor::has_type<first>::value>::observe(
          this_observer, current, previous, r);
      ObserverListImpl<others...>::observe(obs_tuple, current, previous, r);
    }
  };

  template <typename last>
  struct ObserverListImpl<last>
  {
    template <typename T, typename result, typename input>
    static void
    observe(const T&     obs_tuple,
            const input& current,
            const input& previous,
            result&      r)
    {
      const auto& this_observer = std::get<last>(obs_tuple);
      ObserverCaller<result_type_extractor::has_type<last>::value>::observe(
          this_observer, current, previous, r);
    }
  };

  template <>
  struct ObserverListImpl<>
  {
    template <typename T, typename result, typename input>
    static void
    observe(const T&, const input&, result&)
    {
    }
  };
}

template <typename... observers>
struct ObserverList : private detail::Extendable<observers...>
{
private:
  static_assert(not has_duplicates_v<observers...>,
                "same observer type specified several times");

  typedef type_collector_t<result_type_extractor, observers...> results;

  using detail::Extendable<observers...>::tuple;

public:
  template <template <typename...> class R>
  using result_type = unpack_boost_set_as_tparams_t<R, results>;

  using detail::Extendable<observers...>::get;

  template <typename input, typename result_t>
  void
  operator()(const input& current,
             const input& previous,
             result_t&    result) const
  {
    static_assert(all_of_v<observer_traits_checker<observers, input>::value...>,
                  "not all observers support the specified input");

    ObserverListImpl<observers...>::observe(tuple(), current, previous, result);
  }
};

}  // namespace Acts

#endif  // ACTS_OBSERVER_LIST_HPP
