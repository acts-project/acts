// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <tuple>

namespace Acts {

/// Function that allows to zip two ranges to be used in a range-based for loop.
/// This is a very basic implementation for two ranges, but could relatively
/// easily be extended with variadic templates.
/// @tparam RA The first range type
/// @tparam RB The second range type
/// @param ra The first range
/// @param rb The second range
/// @note the behaviour is undefined if the ranges do not have equal range
template <typename RA, typename RB>
auto zip(RA &&ra, RB &&rb) {
  using ItA = decltype(ra.begin());
  using ItB = decltype(rb.begin());

  struct It {
    ItA a;
    ItB b;

    using reference = std::tuple<decltype(*std::declval<ItA>()),
                                 decltype(*std::declval<ItB>())>;

    auto operator++() {
      ++a;
      ++b;
      return *this;
    }

    auto operator!=(const It &other) const { return a != other.a; }

    reference operator*() { return {*a, *b}; }
  };

  struct Zip {
    It b, e;

    auto begin() { return b; }
    auto end() { return e; }
  };

  return Zip{It{ra.begin(), rb.begin()}, It{ra.end(), rb.end()}};
}

}  // namespace Acts
