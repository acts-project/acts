// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <bit>
#include <cassert>

namespace Acts {

namespace Concepts {

/// @brief check types for requirements needed by interpolation
///
/// This helper struct provides compile-time information whether the provided
/// @c Point and @c Value types can be used in the Acts::interpolate function.
///
/// @tparam T type of values to be interpolated
/// @tparam N number of hyper box corners
/// @tparam P1 type for specifying the input point
/// @tparam P2 type for specifying the lower corner of the hyper box
/// @tparam P3 type for specifying the upper corner of the hyper box
template <typename T, std::size_t N, typename P1, typename P2, typename P3>
concept Interpolatable = requires {
  std::has_single_bit(N);
  {
    std::declval<double>() * std::declval<T>() +
        std::declval<double>() * std::declval<T>()
  } -> std::convertible_to<T>;
  { std::declval<P1>()[0] } -> std::convertible_to<double>;
  { std::declval<P2>()[0] } -> std::convertible_to<double>;
  { std::declval<P3>()[0] } -> std::convertible_to<double>;
};

}  // namespace Concepts

/// performs linear interpolation inside a hyper box
///
/// @tparam T type of values to be interpolated
/// @tparam N number of hyper box corners
/// @tparam Point1 type specifying the input point
/// @tparam Point2 type specifying the lower corner of the hyper box
/// @tparam Point3 type specifying the upper corner of the hyper box
///
/// @param point point to which to interpolate
/// @param lowerCorner generalized lower-left corner of hyper box
/// (containing the minima of the hyper box along each dimension)
/// @param upperCorner generalized upper-right corner of hyper box
/// (containing the maxima of the hyper box along each dimension)
/// @param values field values at the hyper box corners sorted in the
/// canonical order defined below.
///
/// @return interpolated value at given point
///
/// @pre @c point must be inside the given hyper box, that
/// is \f$\text{lowerCorner}[i] \le \text{point}[i] \le
/// \text{upperCorner}[i] \quad \forall i=0, \dots, d-1\f$.
///
/// @note
/// - Given @c U and @c V of value type @c T as well as two @c double @c a and
/// @c b, then the following must be a valid expression <tt>a * U + b * V</tt>
/// yielding an object which is (implicitly) convertible to @c T.
/// - All @c Point types must represent d-dimensional points and support
/// coordinate access using @c operator[] which should return a @c double (or a
/// value which is implicitly convertible). Coordinate indices must start at 0.
/// - @c N is the number of hyper box corners which is \f$2^d\f$ where \f$d\f$
/// is the dimensionality of the hyper box. The dimensionality must be
/// consistent with the provided @c Point types.
/// - Definition of the canonical order for sorting the field values: The hyper
/// box corners are numbered according to the following scheme. Each corner is
/// defined by the set of lower/upper boundary limits in each dimension @c i.
/// This can be represented by a binary code (from left to right) where a @c 0
/// stands for a lower bound along this axis and a @c 1 stand for the upper
/// bound along this axis. The left most bit corresponds to the first dimension
/// and the bits to the left correspond to the 2nd, 3rd... dimension. The binary
/// code can be interpreted as integer which gives the number of the
/// corresponding hyper box corner. The field values are ordered according to
/// ascending hyper box corner numbers.<br />
/// As an example assume we have a 3D box with @c lowerCorner = (1,2,3) and @c
/// upperCorner = (4,5,6). The eight corners with their bit patterns and corner
/// numbers are:
///    - (1,2,3): 000 = 0
///    - (1,2,6): 001 = 1
///    - (1,5,3): 010 = 2
///    - (1,5,6): 011 = 3
///    - (4,2,3): 100 = 4
///    - (4,2,6): 101 = 5
///    - (4,5,3): 110 = 6
///    - (4,5,6): 111 = 7
template <typename T, std::size_t N, class Point1, class Point2 = Point1,
          class Point3 = Point2>
T interpolate(const Point1& point, const Point2& lowerCorner,
              const Point3& upperCorner, const std::array<T, N>& values)
  requires Concepts::Interpolatable<T, N, Point1, Point2, Point3>
{
  static constexpr std::size_t D = std::bit_width(N) - 1;
  static constexpr std::size_t I = D - 1;

  // get distance to lower boundary relative to total bin width
  const double f =
      (point[I] - lowerCorner[I]) / (upperCorner[I] - lowerCorner[I]);
  assert(f >= 0 && f <= 1 && "point must be inside the given hyper box");

  std::array<T, N / 2> newValues{};
  for (std::size_t i = 0; i < N / 2; ++i) {
    newValues[i] = (1 - f) * values[2 * i] + f * values[2 * i + 1];
  }

  if constexpr (D == 1) {
    return newValues[0];
  } else {
    return interpolate(point, lowerCorner, upperCorner, newValues);
  }
}

}  // namespace Acts
