// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/detail/interpolation_impl.hpp"

namespace Acts {

/// @brief performs linear interpolation inside a hyper box
///
/// @tparam T     type of values to be interpolated
/// @tparam Point type specifying geometric positions
/// @tparam N     number of hyper box corners
///
/// @param [in] position    position to which to interpolate
/// @param [in] lowerCorner generalized lower-left corner of hyper box
///                         (containing the minima of the hyper box along each
///                         dimension)
/// @param [in] upperCorner generalized upper-right corner of hyper box
///                         (containing the maxima of the hyper box along each
///                         dimension)
/// @param [in] values      field values at the hyper box corners sorted in the
///                         canonical order defined below.
///
/// @return interpolated value at given position
///
/// @pre @c position must describe a position inside the given hyper box, that
///      is \f$\text{lowerCorner}[i] \le \text{position}[i] \le
///      \text{upperCorner}[i] \quad \forall i=0, \dots, d-1\f$.
///
/// @note
/// - Given @c U and @c V of value type @c T as well as two @c double @c a and
/// @c b, then the following must be a valid expression <tt>a * U + b * V</tt>
/// yielding an object which is (implicitly) convertible to @c T.
/// - @c Point must represent a D-dimensional position and support coordinate
/// access using @c operator[]. Coordinate indices must start at 0.
/// - @c N is the number of hyper box corners which is \f$2^D\f$ where \f$D\f$
/// is the dimensionality of the hyper box. The dimensionality must be
/// consistent with the provided @c Point type.
/// - Definition of the canonical order for sorting the field values: The hyper
/// box corners are numbered according to the following scheme. Each corner is
/// defined by the set of lower/upper boundary limits in each dimension @c d.
/// This can be represented by a binary code (from right to left) where a @c 0
/// stands for a lower bound along this axis and a @c 1 stand for the upper
/// bound along this axis. The right most bit corresponds to the first dimension
/// and the bits to the left correspond to the 2nd, 3rd... dimension. The binary
/// code can be interpreted as integer which gives the number of the
/// corresponding hyper box corner. The field values are ordered according to
/// ascending hyper box corner numbers.<br />
/// As an example assume we have a 3D box with @c lowerCorner = (1,2,3) and @c
/// upperCorner = (4,5,6). The eight corners with their bit patterns and corner
/// numbers are:
///    - (1,2,3): 000 = 0
///    - (4,2,3): 001 = 1
///    - (1,5,3): 010 = 2
///    - (4,5,3): 011 = 3
///    - (1,2,6): 100 = 4
///    - (4,2,6): 101 = 5
///    - (1,5,6): 110 = 6
///    - (4,5,6): 111 = 7
template <typename T, class Point, size_t N>
inline T
interpolate(const Point& position,
            const Point& lowerCorner,
            const Point& upperCorner,
            const std::array<T, N>& values)
{
  return detail::interpolate_impl<T, Point, 0u, N>::run(
      position, lowerCorner, upperCorner, values);
}

}  // namespace Acts
