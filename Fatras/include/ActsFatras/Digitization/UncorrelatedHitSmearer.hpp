// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <array>
#include <functional>
#include <utility>

namespace ActsFatras {

/// Smearing function definition for single track parameters.
///
/// The function takes the unsmeared parameter and returns the smeared value and
/// a standard deviation.
///
/// @tparam generator_t The type of the random generator.
template <typename generator_t>
using SingleParameterSmearFunction =
    std::function<Acts::Result<std::pair<double, double>>(double,
                                                          generator_t&)>;

/// Uncorrelated smearing algorithm for fast digitisation of bound parameters.
///
/// @tparam generator_t Random number generator type
/// @tparam kSize Number of smeared parameters
///
/// The smearer takes a single simulated `Hit` and generates a smeared parameter
/// vector and associated covariance matrix.
template <typename generator_t, size_t kSize>
struct BoundParametersSmearer {
  using Scalar = Acts::ActsScalar;
  using ParametersVector = Acts::ActsVector<kSize>;
  using CovarianceMatrix = Acts::ActsSymMatrix<kSize>;
  using Result = Acts::Result<std::pair<ParametersVector, CovarianceMatrix>>;

  /// Parameter indices that will be used to create the smeared measurements.
  std::array<Acts::BoundIndices, kSize> indices{};
  std::array<SingleParameterSmearFunction<generator_t>, kSize> smearFunctions{};

  static constexpr size_t size() { return kSize; }

  /// Generate smeared measured for configured parameters.
  ///
  /// @param rng Random number generator
  /// @param hit Simulated hit
  /// @param surface Local surface on which the hit is smeared
  /// @param geoCtx Geometry context
  /// @retval Smeared parameters vector and associated covariance on success
  /// @retval Error code for failure
  Result operator()(generator_t& rng, const Hit& hit,
                    const Acts::Surface& surface,
                    const Acts::GeometryContext& geoCtx) const {
    // construct full bound parameters. they are probably not all needed, but it
    // is easier to just create them all and then select the requested ones.
    Acts::Result<Acts::BoundVector> boundParamsRes =
        Acts::detail::transformFreeToBoundParameters(hit.position(), hit.time(),
                                                     hit.unitDirection(), 0,
                                                     surface, geoCtx);

    if (!boundParamsRes.ok()) {
      return boundParamsRes.error();
    }

    const auto& boundParams = *boundParamsRes;

    ParametersVector par = ParametersVector::Zero();
    CovarianceMatrix cov = CovarianceMatrix::Zero();
    for (int i = 0; i < static_cast<int>(kSize); ++i) {
      auto res = smearFunctions[i](boundParams[indices[i]], rng);
      if (not res.ok()) {
        return Result::failure(res.error());
      }
      auto [value, stddev] = res.value();
      par[i] = value;
      cov(i, i) = stddev * stddev;
    }

    return Result::success(std::make_pair(par, cov));
  }
};

/// Uncorrelated smearing algorithm for fast digitisation of free parameters.
///
/// @tparam generator_t Random number generator type
/// @tparam kSize Number of smeared parameters
///
/// The smearer takes a single simulated `Hit` and generates a smeared parameter
/// vector and associated covariance matrix.
///
/// @note Uncorrelated smearing of the direction using each components
///   individually is not recommended
template <typename generator_t, size_t kSize>
struct FreeParametersSmearer {
  using Scalar = Acts::ActsScalar;
  using ParametersVector = Acts::ActsVector<kSize>;
  using CovarianceMatrix = Acts::ActsSymMatrix<kSize>;
  using Result = Acts::Result<std::pair<ParametersVector, CovarianceMatrix>>;

  /// Parameter indices that will be used to create the smeared measurements.
  std::array<Acts::FreeIndices, kSize> indices{};
  std::array<SingleParameterSmearFunction<generator_t>, kSize> smearFunctions;

  static constexpr size_t size() { return kSize; }

  /// Generate smeared measured for configured parameters.
  ///
  /// @param rng Random number generator
  /// @param hit Simulated hit
  /// @return Smeared free parameter set wrapped in a Result<...> object
  /// @retval Smeared parameters vector and associated covariance on success
  /// @retval Error code for failure
  Result operator()(generator_t& rng, const Hit& hit) const {
    // construct full free parameters. they are probably not all needed, but it
    // is easier to just create them all and then select the requested ones.
    Acts::FreeVector freeParams;
    freeParams.segment<3>(Acts::eFreePos0) = hit.position();
    freeParams[Acts::eFreeTime] = hit.time();
    freeParams.segment<3>(Acts::eFreeDir0) = hit.unitDirection();
    freeParams[Acts::eFreeQOverP] = 0;

    ParametersVector par = ParametersVector::Zero();
    CovarianceMatrix cov = CovarianceMatrix::Zero();
    for (size_t i = 0; i < kSize; ++i) {
      auto res = smearFunctions[i](freeParams[indices[i]], rng);
      if (not res.ok()) {
        return Result::failure(res.error());
      }
      auto [value, stddev] = res.value();
      par[i] = value;
      cov(i, i) = stddev * stddev;
    }

    return Result::success(std::make_pair(par, cov));
  }
};

}  // namespace ActsFatras
