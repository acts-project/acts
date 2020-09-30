// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <utility>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/Result.hpp"

namespace ActsFatras {

namespace detail {

/// Smearing functions definition:
///
/// @tparam generator_t The type of the random number generator
///
/// - it takes the unsmeared parameter, and a random number generator
/// - it returns the smeared parameter and a covariance
template <typename generator_t>
using SmearFunction = std::function<Acts::Result<std::pair<double, double>>(
    double, generator_t&)>;

/// Single component smearing struct.
///
/// @tparam values_t The type of the values vector
/// @tparam covariance_t The type of the covariance vector
/// @tparam The type of the random number generator
///
/// This is used to unpack the parameter pack smearing
/// @param idx The index of the smearing
/// @param[in,out] values The vector of values to be smeared
/// @param[in,out] covariances The covariance matrix to be set
/// @param sFunction The smearing function to be used
/// @param[in,out] result a recursive smearing result trigger
///        if a single component is not ok, this will be communicated
///        upstream to the caller
template <typename values_t, typename covariance_t, typename generator_t>
struct SingleComponentSmearer {
  void operator()(size_t idx, Eigen::MatrixBase<values_t>& values,
                  Eigen::MatrixBase<covariance_t>& covariances,
                  generator_t& sRandom,
                  const SmearFunction<generator_t>& sFunction,
                  Acts::Result<void>& result) {
    auto sResult = sFunction(values[idx], sRandom);
    if (sResult.ok()) {
      const auto& smeared = sResult.value();
      values[idx] = smeared.first;
      covariances(idx, idx) = smeared.second * smeared.second;
    } else {
      result = Acts::Result<void>(sResult.error());
    }
  }
};

/// Smear the parameters in the given vector with provided smear functions
///
/// @tparam indices_t the Parameter indices type (bound, free)
/// @tparam kParameters Parameter indices pack that defines the used parameter
template <typename indices_t, indices_t... kParameters>
struct ParametersSmearer {
  using StorageSequence = std::make_index_sequence<sizeof...(kParameters)>;

  /// Run the smearing on it
  ///
  /// @tparam values_t The type of the values vector
  /// @tparam covariance_t The type of the covariance vector
  /// @tparam The type of the random number generator
  ///
  /// @param[in,out] values on which smearing will be applied
  /// @param[in,out] covariances which will be filled in place as well
  /// @param[in] sFunctions the smearing functions to be applied
  ///
  /// @return a Result object that may carry a DigitizationError
  template <typename values_t, typename covariance_t, typename generator_t>
  static Acts::Result<void> run(
      Eigen::MatrixBase<values_t>& values,
      Eigen::MatrixBase<covariance_t>& covariances, generator_t& sRandom,
      const std::array<SmearFunction<generator_t>, sizeof...(kParameters)>&
          sFunctions) {
    return runImpl(values, covariances, sRandom, sFunctions, StorageSequence{});
  }

  template <typename values_t, typename covariance_t, typename generator_t,
            std::size_t... kStorage>
  static Acts::Result<void> runImpl(
      Eigen::MatrixBase<values_t>& values,
      Eigen::MatrixBase<covariance_t>& covariances, generator_t& sRandom,
      const std::array<SmearFunction<generator_t>, sizeof...(kParameters)>&
          sFunctions,
      std::index_sequence<kStorage...>) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(values_t, sizeof...(kParameters));
    EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(
        covariance_t, sizeof...(kParameters), sizeof...(kParameters));
    static_assert(sizeof...(kParameters) == sizeof...(kStorage),
                  "Parameters and storage index packs must have the same size");

    Acts::Result<void> result;
    SingleComponentSmearer<values_t, covariance_t, generator_t> scs;
    ((scs(kStorage, values, covariances, sRandom, sFunctions[kStorage],
          result)),
     ...);
    return result;
  }
};

}  // namespace detail
}  // namespace ActsFatras
