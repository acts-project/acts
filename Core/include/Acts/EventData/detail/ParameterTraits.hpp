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

///
/// @brief type for parameters with unrestricted value range
///
struct unbound_parameter {
  static constexpr bool may_modify_value{
      false};  ///< parameter values need no adjustment

  ///
  /// @brief retrieve value for unconstrained parameter value ranges
  ///
  /// @tparam T type of the input parameter
  /// @param input input parameter value
  ///
  /// @return identical input parameter value
  ///
  template <typename T>
  static T getValue(const T& input) {
    return input;
  }

  template <typename T>
  static T getDifference(const T& first, const T& second) {
    return first - second;
  }
};

///
/// @brief type for local parameters bound to a surface
///
struct local_parameter : public unbound_parameter {};

///
/// @brief type for parameter with restricted value range
///
/// This parameter type could be useful to describe parameter with physical
/// meaningful bounds (e.g. radius).
///
/// @tparam T type for boundary value (usually @c double)
/// @tparam MIN pointer to a @c constexpr function returning the lower bound of
/// the value range
/// @tparam MAX pointer to a @c constexpr function returning the upper bound of
/// the value range
///
template <typename T, T (*MIN)(), T (*MAX)()>
struct bound_parameter {
  static constexpr bool may_modify_value{
      true};                      ///< parameter values may need adjustment
  static constexpr T min{MIN()};  ///< lower bound of range
  static constexpr T max{MAX()};  ///< upper bound of range

  ///
  /// @brief retrieve value for constrained parameter value ranges
  ///
  /// @tparam U type of the input parameter
  /// @param input input parameter value
  ///
  /// @return input parameter value cut of at the boundaries @c
  /// bound_parameter<U<MIN<MAX>::min and
  ///         @c bound_parameter<U,MIN,MAX>::max.
  ///
  template <typename U>
  static U getValue(const U& input) {
    return (input > max) ? max : ((input < min) ? min : input);
  }

  template <typename U>
  static U getDifference(const U& first, const U& second) {
    return getValue(first) - getValue(second);
  }
};

///
/// @brief type for parameter with cyclic value range
///
/// This parameter type is useful to e.g. describe angles.
///
/// @tparam T type for boundary value (usually @c double)
/// @tparam MIN pointer to a @c constexpr function returning the lower bound of
/// the value range
/// @tparam MAX pointer to a @c constexpr function returning the upper bound of
/// the value range
///
template <typename T, T (*MIN)(), T (*MAX)()>
struct cyclic_parameter {
  static constexpr bool may_modify_value{
      true};                      ///< parameter values may need adjustment
  static constexpr T min{MIN()};  ///< lower bound of range
  static constexpr T max{MAX()};  ///< upper bound of range

  ///
  /// @brief retrieve value for constrained cyclic parameter value ranges
  ///
  /// @tparam U type of the input parameter
  /// @param input input parameter value
  ///
  /// @return parameter value in the range [@c bound_parameter<U,MIN,MAX>::min,
  ///         @c bound_parameter<U,MIN,MAX>::max] taking into account the cycle
  /// of this
  ///         parameter type.
  ///
  template <typename U>
  static U getValue(const U& input) {
    if (min <= input && input < max) {
      return input;
    } else {
      return input - (max - min) * std::floor((input - min) / (max - min));
    }
  }

  template <typename U>
  static U getDifference(const U& first, const U& second) {
    static constexpr U half_period = (max - min) / 2;
    U tmp = getValue(first) - getValue(second);
    return (tmp < -half_period
                ? tmp + 2 * half_period
                : (tmp > half_period ? tmp - 2 * half_period : tmp));
  }
};

template <BoundIndices>
struct BoundParameterTraits;
template <>
struct BoundParameterTraits<BoundIndices::eBoundLoc0> {
  using type = local_parameter;
};
template <>
struct BoundParameterTraits<BoundIndices::eBoundLoc1> {
  using type = local_parameter;
};
template <>
struct BoundParameterTraits<BoundIndices::eBoundPhi> {
  static constexpr double pMin() { return -M_PI; }
  static constexpr double pMax() { return M_PI; }
  using type = cyclic_parameter<double, pMin, pMax>;
};
template <>
struct BoundParameterTraits<BoundIndices::eBoundTheta> {
  static constexpr double pMin() { return 0; }
  static constexpr double pMax() { return M_PI; }
  using type = bound_parameter<double, pMin, pMax>;
};
template <>
struct BoundParameterTraits<BoundIndices::eBoundQOverP> {
  using type = unbound_parameter;
};
template <>
struct BoundParameterTraits<BoundIndices::eBoundTime> {
  using type = unbound_parameter;
};

template <FreeIndices>
struct FreeParameterTraits;
template <>
struct FreeParameterTraits<FreeIndices::eFreePos0> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreePos1> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreePos2> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreeTime> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreeDir0> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreeDir1> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreeDir2> {
  using type = unbound_parameter;
};
template <>
struct FreeParameterTraits<FreeIndices::eFreeQOverP> {
  using type = unbound_parameter;
};

template <typename indices_t>
struct ParametersSize;
template <>
struct ParametersSize<BoundIndices> {
  static constexpr unsigned int size =
      static_cast<unsigned int>(BoundIndices::eBoundSize);
};
template <>
struct ParametersSize<FreeIndices> {
  static constexpr unsigned int size =
      static_cast<unsigned int>(FreeIndices::eFreeSize);
};

/// Single bound track parameter type for value constrains.
///
/// The singular name is not a typo since this describes individual components.
template <BoundIndices kIndex>
using BoundParameterType = typename detail::BoundParameterTraits<kIndex>::type;

/// Single free track parameter type for value constrains.
///
/// The singular name is not a typo since this describes individual components.
template <FreeIndices kIndex>
using FreeParameterType = typename detail::FreeParameterTraits<kIndex>::type;

/// Access the (Bound/Free)ParameterType through common struct
///
/// @tparam T Parameter indices enum
/// @tparam I Index of @p T
template <typename T, T I>
struct ParameterTypeFor {};

/// Access for @c BoundIndices
template <BoundIndices I>
struct ParameterTypeFor<BoundIndices, I> {
  using type = BoundParameterType<I>;
};

/// Access for @c FreeIndices
template <FreeIndices I>
struct ParameterTypeFor<FreeIndices, I> {
  using type = FreeParameterType<I>;
};

}  // namespace detail
}  // namespace Acts
