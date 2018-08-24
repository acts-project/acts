// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <algorithm>
#include <cmath>

namespace Acts {
///
/// @brief type for parameters with unrestricted value range
///
struct unbound_parameter
{
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
  static T
  getValue(const T& input)
  {
    return input;
  }

  template <typename T>
  static T
  getDifference(const T& first, const T& second)
  {
    return first - second;
  }
};

///
/// @brief type for local parameters bound to a surface
///
struct local_parameter : public unbound_parameter
{
};

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
struct bound_parameter
{
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
  static U
  getValue(const U& input)
  {
    return (input > max) ? max : ((input < min) ? min : input);
  }

  template <typename U>
  static U
  getDifference(const U& first, const U& second)
  {
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
struct cyclic_parameter
{
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
  static U
  getValue(const U& input)
  {
    if (min <= input && input < max) {
      return input;
    } else {
      return input - (max - min) * std::floor((input - min) / (max - min));
    }
  }

  template <typename U>
  static U
  getDifference(const U& first, const U& second)
  {
    static constexpr U half_period = (max - min) / 2;
    U                  tmp         = getValue(first) - getValue(second);
    return (tmp < -half_period
                ? tmp + 2 * half_period
                : (tmp > half_period ? tmp - 2 * half_period : tmp));
  }
};

///
/// @brief parameter traits class
///
/// The user needs to provide a specialization of this class for his parameter
/// definitions.
/// The specialization must include a @c typedef with the name @c
/// parameter_type.
/// This type
/// must fulfill the following requirements:
///   * `parameter_type::getValue` must be a valid expression excepting one
/// parameter value
///     as input and returning a valid parameter value (possibly mapped onto
///     some
/// restricted/
///     cyclic value range).
///   * A static boolean constant with the name `may_modify_value` must exist
/// which indicates
///     whether the `getValue` function may actually return a value different
/// from the output.
///
/// @tparam ParameterPolicy struct or class containing the parameter definitions
/// @tparam parID identifier for the parameter
///
//  template<typename ParameterPolicy,typename ParameterPolicy::par_id_type
//  parID>
//  struct parameter_traits;
}  // namespace Acts