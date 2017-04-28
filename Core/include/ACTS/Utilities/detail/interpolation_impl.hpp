// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>

namespace Acts {

namespace detail {

  /// @brief check types for requirements needed by interpolation
  ///
  /// @tparam Point type for specifying geometric positions
  /// @tparam Value type of values to be interpolated
  ///
  /// This helper struct provides compile-time information whether the provided
  /// @c Point and @c Value types can be used in the Acts::interpolate function.
  ///
  /// The following boolean variable
  /// @code{.cpp}
  /// Acts::detail::can_interpolate<Point,Value>::value
  /// @endcode
  ///
  /// is @c true if both @c Point and @c Value fulfill the type requirements
  /// for being used in the interpolation function, otherwise it is @c false.
  /// This expression can be employed in @c std::enable_if_t to use SFINAE
  /// patterns to enable/disable (member) functions.
  template <typename Point, typename Value>
  struct can_interpolate
  {
    template <typename C>
    static auto
    value_type_test(C* c)
        -> decltype(C(std::declval<double>() * std::declval<C>()
                      + std::declval<double>() * std::declval<C>()),
                    std::true_type());
    template <typename C>
    static std::false_type value_type_test(...);

    template <typename C>
    static auto
    point_type_test(C* c)
        -> decltype(double(std::declval<C>()[0]), std::true_type());
    template <typename C>
    static std::false_type point_type_test(...);

    static const bool value
        = std::is_same<std::true_type,
                       decltype(value_type_test<Value>(nullptr))>::value
        and std::is_same<std::true_type,
                         decltype(point_type_test<Point>(nullptr))>::value;
  };

  /// @brief helper struct for performing multi-dimensional linear interpolation
  ///
  /// @tparam T     type of values to be interpolated
  /// @tparam Point type specifying geometric positions
  /// @tparam D     current dimension on which to perform reduction
  /// @tparam N     number of hyper box corners
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
  template <typename T, class Point, size_t D, size_t N>
  struct interpolate_impl;

  /// @cond
  // recursive implementation of linear interpolation in multiple dimensions
  template <typename T, class Point, size_t D, size_t N>
  struct interpolate_impl
  {
    static T
    run(const Point& pos,
        const Point& lowerLeft,
        const Point& upperRight,
        const std::array<T, N>& fields)
    {
      // get distance to lower boundary relative to total bin width
      const double f = (pos[D] - lowerLeft[D]) / (upperRight[D] - lowerLeft[D]);

      std::array<T, (N >> 1)> newFields;
      for (size_t i     = 0; i < N / 2; ++i)
        newFields.at(i) = (1 - f) * fields.at(2 * i) + f * fields.at(2 * i + 1);

      return interpolate_impl<T, Point, D + 1, (N >> 1)>::run(
          pos, lowerLeft, upperRight, newFields);
    }
  };

  // simple linear interpolation in 1D
  template <typename T, class Point, size_t D>
  struct interpolate_impl<T, Point, D, 2u>
  {
    static T
    run(const Point& pos,
        const Point& lowerLeft,
        const Point& upperRight,
        const std::array<T, 2u>& fields)
    {
      // get distance to lower boundary relative to total bin width
      const double f = (pos[D] - lowerLeft[D]) / (upperRight[D] - lowerLeft[D]);

      return (1 - f) * fields.at(0) + f * fields.at(1);
    }
  };
  /// @endcond
}  // namespace detail

}  // namespace Acts
