// This file is part of the Acts project.
//
// Copyright (C) 2016 Acts project team
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
      typedef typename bm::set<args...>::type type;
    };

    template <typename T, typename R>
    struct fold2tuple;

    template <typename... args, typename next>
    struct fold2tuple<std::tuple<args...>, next>
    {
      typedef std::tuple<next, args...> type;
    };

    template <typename S>
    struct boost_set2tuple
    {
      typedef typename bm::fold<S, std::tuple<>, fold2tuple<bp::_1, bp::_2>>::type type;
    };

    template <typename s, typename v>
    struct boost_set_merger
    {
      typedef typename bm::fold<v, s, bm::insert<bp::_1, bp::_2>>::type                       unique_types;
      typedef typename bm::fold<unique_types, std::tuple<>, fold2tuple<bp::_1, bp::_2>>::type flat_types;
      typedef typename tuple2boost_set<flat_types>::type                                      type;
    };

    template <template <typename...> class R, typename S>
    struct boost_set_as_tparams
    {
      typedef typename boost_set2tuple<S>::type                as_tuple;
      typedef typename boost_set_as_tparams<R, as_tuple>::type type;
    };

    template <template <typename...> class R, typename... args>
    struct boost_set_as_tparams<R, std::tuple<args...>>
    {
      typedef R<args...> type;
    };
    // clang-format on
  }  // end of anonymous namespace

  template <typename s, typename v>
  using boost_set_merger_t = typename boost_set_merger<s, v>::type;

  template <template <typename...> class R, typename S>
  using boost_set_as_tparams_t = typename boost_set_as_tparams<R, S>::type;

}  // namespace detail

}  // namespace Acts