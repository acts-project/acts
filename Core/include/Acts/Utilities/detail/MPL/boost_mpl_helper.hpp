// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <boost/mpl/fold.hpp>
#include <boost/mpl/inserter.hpp>
#include <boost/mpl/set.hpp>
#include <tuple>

namespace Acts {

namespace detail {

  namespace {
    // clang-format off
    namespace bm = boost::mpl;
    namespace bp = boost::mpl::placeholders;

    template <typename seq>
    struct tuple2boost_set;

    template <typename... args>
    struct tuple2boost_set<std::tuple<args...>>
    {
      using type = typename bm::set<args...>::type;
    };

    template <typename T, typename R>
    struct fold2tuple;

    template <typename... args, typename next>
    struct fold2tuple<std::tuple<args...>, next>
    {
      using type = std::tuple<next, args...>;
    };

    template <typename S>
    struct boost_set2tuple
    {
      using type = typename bm::fold<S, std::tuple<>, fold2tuple<bp::_1, bp::_2>>::type;
    };

    template <typename s, typename v>
    struct boost_set_merger
    {
      using unique_types = typename bm::fold<v, s, bm::insert<bp::_1, bp::_2>>::type;
      using flat_types = typename bm::fold<unique_types, std::tuple<>, fold2tuple<bp::_1, bp::_2>>::type;
      using type = typename tuple2boost_set<flat_types>::type;
    };

    template <template <typename...> class R, typename S>
    struct boost_set_as_tparams
    {
      using as_tuple = typename boost_set2tuple<S>::type;
      using type = typename boost_set_as_tparams<R, as_tuple>::type;
    };

    template <template <typename...> class R, typename... args>
    struct boost_set_as_tparams<R, std::tuple<args...>>
    {
      using type = R<args...>;
    };
    // clang-format on
  }  // end of anonymous namespace

  template <typename s, typename v>
  using boost_set_merger_t = typename boost_set_merger<s, v>::type;

  template <template <typename...> class R, typename S>
  using boost_set_as_tparams_t = typename boost_set_as_tparams<R, S>::type;

}  // namespace detail

}  // namespace Acts