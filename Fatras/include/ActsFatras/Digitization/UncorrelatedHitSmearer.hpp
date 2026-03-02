// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"
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
template <typename generator_t, std::size_t kSize>
struct BoundParametersSmearer {
  /// Type alias for parameter vector of dimension kSize
  using ParametersVector = Acts::Vector<kSize>;
  /// Type alias for covariance matrix of dimension kSize x kSize
  using CovarianceMatrix = Acts::SquareMatrix<kSize>;
  /// Type alias for smearing result containing parameters and covariance
  using Result = Acts::Result<std::pair<ParametersVector, CovarianceMatrix>>;

  /// Parameter indices that will be used to create the smeared measurements.
  std::array<Acts::BoundIndices, kSize> indices{};
  /// Array of smearing functions for each measurement parameter
  std::array<SingleParameterSmearFunction<generator_t>, kSize> smearFunctions{};
  /// Array of flags to force positive values after smearing
  std::array<bool, kSize> forcePositive = {};
  /// Maximum number of retries for generating positive values when forced
  std::size_t maxRetries = 0;

  /// Get the number of bound parameters that will be smeared
  /// @return Number of bound parameters (kSize)
  static constexpr std::size_t size() { return kSize; }

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
    // We use the thickness of the detector element as tolerance, because Geant4
    // treats the Surfaces as volumes and thus it is not ensured, that each hit
    // lies exactly on the Acts::Surface
    const auto tolerance = surface.isSensitive() ? surface.thickness()
                                                 : Acts::s_onSurfaceTolerance;

    // construct full bound parameters. they are probably not all needed, but it
    // is easier to just create them all and then select the requested ones.
    Acts::Result<Acts::BoundVector> boundParamsRes =
        Acts::transformFreeToBoundParameters(hit.position(), hit.time(),
                                             hit.direction(), 0, surface,
                                             geoCtx, tolerance);

    if (!boundParamsRes.ok()) {
      return boundParamsRes.error();
    }

    const auto& boundParams = *boundParamsRes;
    Acts::BoundVector smearedBoundParams = boundParams;

    for (std::size_t k = 0; k < maxRetries + 1; ++k) {
      ParametersVector par = ParametersVector::Zero();
      CovarianceMatrix cov = CovarianceMatrix::Zero();
      for (std::size_t i = 0; i < kSize; ++i) {
        auto res = smearFunctions[i](boundParams[indices[i]], rng);
        if (!res.ok()) {
          return Result::failure(res.error());
        }
        auto [value, stddev] = res.value();
        par[i] = value;
        if (forcePositive[i]) {
          par[i] = std::abs(value);
        }
        smearedBoundParams[indices[i]] = par[i];
        cov(i, i) = stddev * stddev;
      }

      if (!surface.insideBounds(smearedBoundParams.head<2>())) {
        continue;
      }

      return Result::success(std::make_pair(par, cov));
    }

    return Result::failure(DigitizationError::MaximumRetriesExceeded);
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
template <typename generator_t, std::size_t kSize>
struct FreeParametersSmearer {
  /// Type alias for parameter vector of dimension kSize
  using ParametersVector = Acts::Vector<kSize>;
  /// Type alias for covariance matrix of dimension kSize x kSize
  using CovarianceMatrix = Acts::SquareMatrix<kSize>;
  /// Type alias for smearing result containing parameters and covariance
  using Result = Acts::Result<std::pair<ParametersVector, CovarianceMatrix>>;

  /// Parameter indices that will be used to create the smeared measurements.
  std::array<Acts::FreeIndices, kSize> indices{};
  /// Array of smearing functions for each free parameter
  std::array<SingleParameterSmearFunction<generator_t>, kSize> smearFunctions;

  /// Get the number of free parameters that will be smeared
  /// @return Number of free parameters (kSize)
  static constexpr std::size_t size() { return kSize; }

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
    freeParams.segment<3>(Acts::eFreeDir0) = hit.direction();
    freeParams[Acts::eFreeQOverP] = 0;

    ParametersVector par = ParametersVector::Zero();
    CovarianceMatrix cov = CovarianceMatrix::Zero();
    for (std::size_t i = 0; i < kSize; ++i) {
      auto res = smearFunctions[i](freeParams[indices[i]], rng);
      if (!res.ok()) {
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
