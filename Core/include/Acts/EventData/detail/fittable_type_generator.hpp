// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// boost include(s)
//#include <boost/mpl/vector.hpp>
//#include <boost/variant.hpp>
#include "Acts/Utilities/ParameterDefinitions.hpp"

#include <type_traits>

#include <boost/hana/append.hpp>
#include <boost/hana/integral_constant.hpp>
#include <boost/hana/range.hpp>
#include <boost/hana/remove_at.hpp>
#include <boost/hana/transform.hpp>
#include <boost/hana/tuple.hpp>

namespace Acts {
// forward declaration
template <typename identifier_t, ParID_t... params>
class Measurement;

/// @cond detail
namespace detail {

  namespace hana = boost::hana;

  template <size_t W>
  constexpr auto
  unique_ordered_sublists()
  {
    using namespace hana::literals;
    constexpr auto combinations = hana::make_tuple(hana::make_tuple());
    constexpr auto w_range
        = hana::to_tuple(hana::make_range(0_c, hana::size_c<W>));
    constexpr auto comb2
        = hana::fold_left(w_range, combinations, [](auto state, auto i) {
            auto mapped = hana::transform(
                state, [i](auto c_i) { return hana::append(c_i, i); });
            return hana::concat(state, mapped);
          });
    return hana::remove_at(comb2, 0_c);
  }

  template <template <ParID_t...> class meas_meta, size_t W>
  constexpr auto
  type_generator()
  {
    constexpr auto sublists       = unique_ordered_sublists<W>();
    constexpr auto measurements_h = hana::transform(sublists, [](auto s) {
      return hana::unpack(s, [](auto... i) {
        return hana::type_c<typename meas_meta<ParID_t(decltype(i)::value)...>::type>;
      });
    });
    return measurements_h;
  }

  template <template <ParID_t...> class meas_meta, size_t W>
  using type_generator_t =
      typename decltype(hana::unpack(type_generator<meas_meta, W>(),
                                     hana::template_<std::variant>))::type;

  /// @endcond
}  // namespace detail
/// @endcond
}  // namespace Acts
