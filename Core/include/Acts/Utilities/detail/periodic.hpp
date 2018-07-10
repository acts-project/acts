// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cmath>

namespace Acts {
namespace detail {

  /// Wrap a periodic value back into the nominal range.
  template <typename T>
  inline T
  wrap_periodic(T value, T start, T range)
  {
    using std::floor;
    // only wrap if really necessary
    T diff = value - start;
    return ((0 <= diff) && (diff < range))
        ? value
        : (value - range * floor(diff / range));
  }

  /// Calculate the equivalent angle in the [0, 2*pi) range.
  template <typename T>
  inline T
  radian_pos(T x)
  {
    return wrap_periodic<T>(x, T(0), T(2 * M_PI));
  }

  /// Calculate the equivalent angle in the [-pi, pi) range.
  template <typename T>
  inline T
  radian_sym(T x)
  {
    return wrap_periodic<T>(x, T(-M_PI), T(2 * M_PI));
  }

}  // namespace detail
}  // namespace Acts