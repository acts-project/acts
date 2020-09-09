// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/ParameterDefinitions.hpp"

#include <algorithm>
#include <cmath>

namespace Acts {
namespace detail {

/// Traits class for an unrestricted parameter.
struct UnrestrictedParameterTraits {
  /// Values need no adjustment.
  static constexpr bool may_modify_value = false;

  /// Get the corrected value within the limits. This is a no-op here.
  template <typename value_t>
  static constexpr const value_t& getValue(const value_t& value) {
    return value;
  }
  /// Compute the difference between two values.
  template <typename value_t>
  static constexpr value_t getDifference(const value_t& lhs,
                                         const value_t& rhs) {
    return lhs - rhs;
  }
};

/// Traits class for a parameter with a restricted value range.
///
/// @tparam limits_t a type with static `lowest()` and `max()` member functions
///
/// This parameter type could be useful to describe parameter with physical
/// meaningful bounds (e.g. radius).
template <typename limits_t>
struct RestrictedParameterTraits {
  /// Parameter values may need adjustment.
  static constexpr bool may_modify_value = true;
  /// Lower bound of range.
  static constexpr double min = limits_t::lowest();
  /// Upper bound of range.
  static constexpr double max = limits_t::max();

  /// Get the corrected value within the limits.
  template <typename value_t>
  static constexpr value_t getValue(const value_t& value) {
    return std::clamp(value, static_cast<value_t>(min),
                      static_cast<value_t>(max));
  }
  /// Compute the difference between two values with limit handling.
  template <typename value_t>
  static constexpr value_t getDifference(const value_t& lhs,
                                         const value_t& rhs) {
    return getValue(lhs) - getValue(rhs);
  }
};

/// Traits class for a parameter with a cyclic value range.
///
/// @tparam limits_t a type with static `lowest()` and `max()` member functions
///
/// This parameter type is useful to e.g. describe angles.
template <typename limits_t>
struct CyclicParameterTraits {
  /// Parameter values may need adjustment.
  static constexpr bool may_modify_value = true;
  /// Lower bound of range.
  static constexpr double min = limits_t::lowest();
  /// Upper bound of range.
  static constexpr double max = limits_t::max();

  /// Get the corrected value folded into the central range.
  template <typename value_t>
  static constexpr value_t getValue(const value_t& value) {
    if ((min <= value) and (value < max)) {
      return value;
    } else {
      return value - (max - min) * std::floor((value - min) / (max - min));
    }
  }
  /// Compute the smallest equivalent difference when using the periodicity.
  template <typename value_t>
  static constexpr value_t getDifference(const value_t& lhs,
                                         const value_t& rhs) {
    constexpr value_t range = (max - min);
    value_t delta = getValue(lhs) - getValue(rhs);
    if ((2 * delta) < -range) {
      return delta + range;
    } else if (range < (2 * delta)) {
      return delta - range;
    } else {
      return delta;
    }
  }
};

// Traits types for bound parameters.
//
// These types should not be used directly but only via `ParameterTraits` below.

struct PhiBoundParameterLimits {
  // use function names consistent w/ std::numeric_limits
  static constexpr double lowest() { return -M_PI; }
  static constexpr double max() { return M_PI; }
};
struct ThetaBoundParameterLimits {
  // use function names consistent w/ std::numeric_limits
  static constexpr double lowest() { return 0; }
  static constexpr double max() { return M_PI; }
};

// all parameters not explicitely specified are unrestricted
template <BoundIndices>
struct BoundParameterTraits : public UnrestrictedParameterTraits {};
template <>
struct BoundParameterTraits<BoundIndices::eBoundPhi>
    : public CyclicParameterTraits<PhiBoundParameterLimits> {};
template <>
struct BoundParameterTraits<BoundIndices::eBoundTheta>
    : public RestrictedParameterTraits<ThetaBoundParameterLimits> {};

// The separate traits implementation structs are needed to generate a
// compile-time error in case a non-supported enum type/ index combinations is
// used.
template <typename index_t, index_t kIndex>
struct ParameterTraitsImpl;
template <BoundIndices kIndex>
struct ParameterTraitsImpl<BoundIndices, kIndex> {
  using Type = BoundParameterTraits<kIndex>;
};
template <FreeIndices kIndex>
struct ParameterTraitsImpl<FreeIndices, kIndex> {
  // all free parameters components are unrestricted
  using Type = UnrestrictedParameterTraits;
};

/// Parameter traits for an index from one of the indices enums.
///
/// @tparam index_t Parameter indices enum
/// @tparam kIndex A specific parameter index
///
/// This type resolves directly to one of the parameter traits classes defined
/// above and allows for uniform access.
template <typename index_t, index_t kIndex>
using ParameterTraits = typename ParameterTraitsImpl<index_t, kIndex>::Type;

// The separate size implementation structs are needed to generate a
// compile-time error in case a non-supported enum type combinations is used.
// Template variables can not have an undefined default case.
template <typename indices_t>
struct ParametersSizeImpl;
template <>
struct ParametersSizeImpl<BoundIndices> {
  static constexpr unsigned int kSize =
      static_cast<unsigned int>(BoundIndices::eBoundSize);
};
template <>
struct ParametersSizeImpl<FreeIndices> {
  static constexpr unsigned int kSize =
      static_cast<unsigned int>(FreeIndices::eFreeSize);
};

/// The maximum parameters vector size definable by an index enum.
template <typename indices_t>
constexpr unsigned int kParametersSize = ParametersSizeImpl<indices_t>::kSize;

}  // namespace detail
}  // namespace Acts
