// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
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

/// Ensure both phi and theta direction angles are within the allowed range.
///
/// @param[in] phi Transverse direction angle
/// @param[in] theta Longitudinal direction angle
/// @return pair<phi,theta> containing the updated angles
///
/// The phi angle is truly cyclic, i.e. all values outside the nominal range
/// [-pi,pi) have a corresponding value inside nominal range, independent from
/// the theta angle. The theta angle is more complicated. Imagine that the two
/// angles describe a position on the unit sphere. If theta moves outside its
/// nominal range [0,pi], we are moving over one of the two poles of the unit
/// sphere along the great circle defined by phi. The angles still describe a
/// valid position on the unit sphere, but to describe it with angles within
/// their nominal range, both phi and theta need to be updated; when moving over
/// the poles, phi needs to be flipped by 180degree to allow theta to remain
/// within its nominal range.
template <typename T>
inline std::pair<T, T> normalizePhiTheta(T phi, T theta) {
  // wrap to [0,2pi). while the nominal range of theta is [0,pi], it is
  // periodic, i.e. describes identical positions, in the full [0,2pi) range.
  // moving it first to the periodic range simplifies further steps as the
  // possible range of theta becomes fixed.
  theta = radian_pos(theta);
  if (M_PI < theta) {
    // theta is in the second half of the great circle and outside its nominal
    // range. need to change both phi and theta to be within range.
    phi += M_PI;
    theta = 2 * M_PI - theta;
  }
  return {radian_sym(phi), theta};
}

}  // namespace detail
}  // namespace Acts
