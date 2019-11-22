// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
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
inline T wrap_periodic(T value, T start, T range) {
  using std::floor;
  // only wrap if really necessary
  T diff = value - start;
  return ((0 <= diff) && (diff < range))
             ? value
             : (value - range * floor(diff / range));
}

/// Calculate the equivalent angle in the [0, 2*pi) range.
template <typename T>
inline T radian_pos(T x) {
  return wrap_periodic<T>(x, T(0), T(2 * M_PI));
}

/// Calculate the equivalent angle in the [-pi, pi) range.
template <typename T>
inline T radian_sym(T x) {
  return wrap_periodic<T>(x, T(-M_PI), T(2 * M_PI));
}

/// Calculates the equivalent angles phi and theta
/// in the [-pi, pi) and [0, pi] ranges by ensuring
/// the correct theta bounds
template <typename T>
inline std::pair<T, T> ensureThetaBounds(T phi, T theta) {
  T tmpPhi = radian_sym(phi);

  T tmpTht = std::fmod(theta, 2 * M_PI);
  if (tmpTht < -M_PI) {
    tmpTht = std::abs(tmpTht + 2 * M_PI);
  } else if (tmpTht < 0) {
    tmpTht *= -1;
    tmpPhi += M_PI;
    tmpPhi = tmpPhi > M_PI ? tmpPhi - 2 * M_PI : tmpPhi;
  }
  if (tmpTht > M_PI) {
    tmpTht = 2 * M_PI - tmpTht;
    tmpPhi += M_PI;
    tmpPhi = tmpPhi > M_PI ? (tmpPhi - 2 * M_PI) : tmpPhi;
  }

  return std::pair<T, T>(tmpPhi, tmpTht);
}

}  // namespace detail
}  // namespace Acts