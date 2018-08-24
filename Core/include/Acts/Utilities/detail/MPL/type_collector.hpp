// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <boost/mpl/set.hpp>
#include <type_traits>
#include "Acts/Utilities/detail/MPL/boost_mpl_helper.hpp"
namespace bm = boost::mpl;

#define HAS_TYPE(x)                                                            \
  template <typename T, typename = void>                                       \
  struct has_type : public std::false_type                                     \
  {                                                                            \
  };                                                                           \
                                                                               \
  template <typename T>                                                        \
  struct has_type<T,                                                           \
                  typename std::enable_if<(sizeof(typename T::x) > 0),         \
                                          void>::type> : public std::true_type \
  {                                                                            \
  };

namespace Acts {

namespace detail {

  namespace {

    struct result_type_extractor
    {
    private:
      template <typename T>
      struct extractor
      {
        using type = typename T::result_type;
      };

    public:
      HAS_TYPE(result_type)

      template <typename T>
      using type = typename extractor<T>::type;
    };

    struct action_type_extractor
    {
    private:
      template <typename T>
      struct extractor
      {
        using type = typename T::action_type;
      };

    public:
      HAS_TYPE(action_type)

      template <typename T>
      using type = typename extractor<T>::type;
    };

    template <typename sequence, typename ex, typename T, bool has_type = false>
    struct type_inserter
    {
      using type = sequence;
    };

    template <typename sequence, typename ex, typename T>
    struct type_inserter<sequence, ex, T, true>
    {
      using type =
          typename bm::insert<sequence, typename ex::template type<T>>::type;
    };

    template <typename sequence, typename ex, typename... traits>
    struct type_collector_impl;

    template <typename sequence, typename ex>
    struct type_collector_impl<sequence, ex>
    {
      using type = sequence;
    };

    template <typename sequence,
              typename ex,
              typename first,
              typename... others>
    struct type_collector_impl<sequence, ex, first, others...>
    {
      using new_seq =
          typename type_inserter<sequence,
                                 ex,
                                 first,
                                 ex::template has_type<first>::value>::type;
      using type = typename type_collector_impl<new_seq, ex, others...>::type;
    };

    template <typename sequence, typename ex, typename last>
    struct type_collector_impl<sequence, ex, last>
    {
      using type =
          typename type_inserter<sequence,
                                 ex,
                                 last,
                                 ex::template has_type<last>::value>::type;
    };

    template <typename extractor, typename... traits>
    struct type_collector
    {
      using found =
          typename type_collector_impl<bm::set<>, extractor, traits...>::type;
      using type = boost_set_merger_t<found, bm::set<>>;
    };
  }  // end of anonymous namespace

  template <typename T>
  constexpr bool has_result_type_v = result_type_extractor::has_type<T>::value;

  template <typename T>
  using result_type_t = typename result_type_extractor::type<T>;

  template <typename T>
  constexpr bool has_action_type_v = action_type_extractor::has_type<T>::value;

  template <typename T>
  using action_type_t = typename action_type_extractor::type<T>;

  template <typename extractor, typename... traits>
  using type_collector_t = typename type_collector<extractor, traits...>::type;
}  // namespace detail

}  // namespace Acts
