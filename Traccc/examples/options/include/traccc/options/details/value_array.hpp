// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// System include(s).
#include <array>
#include <iosfwd>

namespace traccc::opts {

/// A fixed number of values as one user option.
///
/// @note Implemented as a subclass so it is distinct from `std::array`
///   and we can provide overloads in the same namespace.
///
template <typename TYPE, std::size_t SIZE>
class value_array : public std::array<TYPE, SIZE> {
 public:
  /// Helper function for normalizing the values in an array
  ///
  /// @param factor The factor to normalize by
  /// @return The normalized array
  ///
  value_array& operator*=(TYPE factor);
  /// Helper function for normalizing the values in an array
  ///
  /// @param factor The factor to normalize by
  /// @return The normalized array
  ///
  value_array& operator/=(TYPE factor);

};  // class value_array

/// Extract a fixed number of values from an input of the form 'x:y:z'.
///
/// @note If the values would be separated by whitespace, negative values
///   and additional command line both start with `-` and would be
///   undistinguishable.
///
/// @param is The input stream to read from.
/// @param values The array to fill with the values.
/// @return The input stream for chaining.
///
template <typename TYPE, std::size_t SIZE>
std::istream& operator>>(std::istream& is, value_array<TYPE, SIZE>& values);

/// Print a fixed number of values as `x:y:z`.
///
/// @param os The output stream to write to.
/// @param values The array to print.
/// @return The output stream for chaining.
///
template <typename TYPE, std::size_t SIZE>
std::ostream& operator<<(std::ostream& os,
                         const value_array<TYPE, SIZE>& values);

/// Helper function for normalizing the values in an array
///
/// @param values The array to normalize
/// @param factor The factor to normalize by
/// @return The normalized array
///
template <typename TYPE, std::size_t SIZE>
value_array<TYPE, SIZE> operator*(const value_array<TYPE, SIZE>& values,
                                  TYPE factor);

/// Helper function for normalizing the values in an array
///
/// @param values The array to normalize
/// @param factor The factor to normalize by
/// @return The normalized array
///
template <typename TYPE, std::size_t SIZE>
value_array<TYPE, SIZE> operator/(const value_array<TYPE, SIZE>& values,
                                  TYPE factor);

}  // namespace traccc::opts

// Include the implementation.
#include "value_array.ipp"
