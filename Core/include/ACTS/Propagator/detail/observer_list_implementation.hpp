#ifndef ACTS_OBSERVER_LIST_IMPLEMENTATION_HPP
#define ACTS_OBSERVER_LIST_IMPLEMENTATION_HPP 1

#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  namespace {
    template <bool has_result = true>
    struct observer_caller
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
            r.template get<detail::result_type_t<observer>>());
      }
    };

    template <>
    struct observer_caller<false>
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
  }  // end of anonymous namespace

  template <typename... observers>
  struct observer_list_impl;

  template <typename first, typename... others>
  struct observer_list_impl<first, others...>
  {
    template <typename T, typename result, typename input>
    static void
    observe(const T&     obs_tuple,
            const input& current,
            const input& previous,
            result&      r)
    {
      constexpr bool has_result    = has_result_type_v<first>;
      const auto&    this_observer = std::get<first>(obs_tuple);
      observer_caller<has_result>::observe(this_observer, current, previous, r);
      observer_list_impl<others...>::observe(obs_tuple, current, previous, r);
    }
  };

  template <typename last>
  struct observer_list_impl<last>
  {
    template <typename T, typename result, typename input>
    static void
    observe(const T&     obs_tuple,
            const input& current,
            const input& previous,
            result&      r)
    {
      constexpr bool has_result    = has_result_type_v<last>;
      const auto&    this_observer = std::get<last>(obs_tuple);
      observer_caller<has_result>::observe(this_observer, current, previous, r);
    }
  };

  template <>
  struct observer_list_impl<>
  {
    template <typename T, typename result, typename input>
    static void
    observe(const T&, const input&, const input&, result&)
    {
    }
  };
}  // namespace detail

}  // namespace Acts
#endif  // ACTS_OBSERVER_LIST_IMPLEMENTATION_HPP
