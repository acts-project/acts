// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <exception>
#include <system_error>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/Result.hpp"

namespace ActsFatras {

/// Smearing functions definition, it retuns the smeared parameter and the
/// (approximated) error for the covariance entry
using SmearFunction = std::function<Acts::Result<std::pair<double, double>>()>;

/// Parameter smearer for fast digitisation that can produce surface bound
/// and free measurements.
///
/// This smearer only supports uncorrelated smearing of parameters
class UncorrelatedHitSmearer {
 public:
  /// Shorthand for Measurement on Surfaces (bound parameterisation)
  template <typename hit_t, Acts::BoundIndices... params>
  using SurfaceMeasurement =
      Acts::Measurement<hit_t, Acts::BoundIndices, params...>;

  /// Shorthand for Measurement in Volumes (free parameterisation)
  template <typename hit_t, Acts::FreeIndices... params>
  using VolumeMeasurement =
      Acts::Measurement<hit_t, Acts::FreeIndices, params...>;

  /// Generic implementation of a smearing meathod
  /// for surface based Measurements. This can be applied with custom
  /// smearing functions and to variant type measurements.
  ///
  /// @tparam hit_t definition requires position(), unitDirection(), time()
  /// @tparam params parameter pack describing the parameters to smear
  ///
  /// @param gctx The Geometry context (alignment, etc.)
  /// @param hit The Fatras hit to be smeared (will be source linked)
  /// @param sf The Surface to which this hit is attached
  /// @param sFunctions The smearing functions that are applied
  ///
  /// @return  A surface based Measurement, or a null result
  template <typename hit_t, Acts::BoundIndices... params>
  Acts::Result<SurfaceMeasurement<hit_t, params...>> createSurfaceMeasurement(
      const Acts::GeometryContext& gctx, const hit_t& hit,
      const Acts::Surface& sf,
      const std::array<SmearFunction, sizeof...(params)>& sFunctions) const {
    using Result = Acts::Result<SurfaceMeasurement<hit_t, params...>>;

    auto dir = hit.unitDirection();
    auto gltResult = sf.globalToLocal(gctx, hit.position(), dir);
    if (not gltResult.ok()) {
      throw Result(Acts::SurfaceError::GlobalPositionNotOnSurface);
    }
    const auto& localHit = gltResult.value();

    Acts::BoundVector boundVector;
    boundVector[Acts::eBoundLoc0] = localHit[Acts::eBoundLoc0];
    boundVector[Acts::eBoundLoc1] = localHit[Acts::eBoundLoc1];
    boundVector[Acts::eBoundPhi] = Acts::VectorHelpers::phi(dir);
    boundVector[Acts::eBoundTheta] = Acts::VectorHelpers::theta(dir);
    boundVector[Acts::eBoundTime] = hit.time();

    typename SurfaceMeasurement<hit_t, params...>::ParameterVector vector;
    vector.setZero();

    typename SurfaceMeasurement<hit_t, params...>::CovarianceMatrix covariance;
    covariance.setZero();

    auto sResult =
        smear(boundVector, vector, covariance, sFunctions, params...);
    if (not sResult.ok()) {
      return Result(std::error_code());
    }

    return Result(SurfaceMeasurement<hit_t, params...>(
        sf.getSharedPtr(), hit, std::move(covariance), vector));
  }

  /// Generic implementation of a smearing meathod
  /// for volume based Measurements. This can be applied with custom
  /// smearing functions and to variant type measurements.
  ///
  /// @tparam hit_t definition requires position(), unitDirection(), time()
  /// @tparam params parameter pack describing the parameters to smear
  ///
  /// @param hit The Fatras hit to be smeared (will be source linked)
  /// @param volume The Volume in which this measurement was measured
  /// @param sFunctions The smearing functions that are applied
  ///
  /// @return  A surface based Measurement, or a null result
  template <typename hit_t, Acts::FreeIndices... params>
  Acts::Result<VolumeMeasurement<hit_t, params...>> createVolumeMeasurement(
      const hit_t& hit, std::shared_ptr<const Acts::Volume> volume,
      const std::array<SmearFunction, sizeof...(params)>& sFunctions) const {
    using Result = Acts::Result<VolumeMeasurement<hit_t, params...>>;

    auto pos = hit.position();
    auto dir = hit.unitDirection();

    Acts::FreeVector freeVector;
    freeVector[Acts::eFreePos0] = pos.x();
    freeVector[Acts::eFreePos1] = pos.y();
    freeVector[Acts::eFreePos2] = pos.z();
    freeVector[Acts::eFreeTime] = hit.time();
    freeVector[Acts::eFreeDir0] = dir.x();
    freeVector[Acts::eFreeDir1] = dir.y();
    freeVector[Acts::eFreeDir2] = dir.z();

    typename VolumeMeasurement<hit_t, params...>::ParameterVector vector;
    vector.setZero();

    typename VolumeMeasurement<hit_t, params...>::CovarianceMatrix covariance;
    covariance.setZero();

    auto sResult = smear(freeVector, vector, covariance, sFunctions, params...);
    if (not sResult.ok()) {
      return Result(sResult.error());
    }

    return Result(VolumeMeasurement<hit_t, params...>(
        volume, hit, std::move(covariance), vector));
  }

 private:
  /// Apply the recursive smearing on the parameter pack,
  /// the base local position or time is used as unsmeared
  ///
  /// @tparam full_vec_t The full parameter vector, unsmeared
  /// @tparam vec_t The parameter vector, static size known at compile time
  /// @tparam cov_t The covariance materix, static size known at compile time
  /// @tparam smear_t The smearing funtion array, static size known at c.t.
  /// @tparam param_t The first of the pack
  /// @tparam paramts_t The remaining pack
  ///
  /// @param fVector The original unsmeared hit information (full dimensional)
  /// @param vector The to be filled parameters vector
  /// @param covariance The to be filled covariance matrix
  /// @param smearing The array of smearing functions
  /// @param first The first of the pack
  /// @param rest The remaining pack
  ///
  /// @return boolean to indiciate it was successful
  template <typename full_vec_t, typename vec_t, typename cov_t,
            typename smear_t, typename param_t, typename... params_t>
  Acts::Result<void> smear(const full_vec_t& fVector, vec_t& vector,
                           cov_t& covariance, smear_t& smearing,
                           const param_t& first,
                           params_t const&... rest) const {
    const size_t par = vector.rows() - sizeof...(rest) - 1;
    auto sResult = smearing[par]();
    if (sResult.ok()) {
      const auto& sr = sResult.value();
      vector[par] = fVector[first] + sr.first;
      covariance(par, par) = sr.second * sr.second;

      if constexpr (sizeof...(rest) > 0) {
        if (not smear(fVector, vector, covariance, smearing, rest...).ok()) {
          return Acts::Result<void>::failure(sResult.error());
        }
      }
      return Acts::Result<void>::success();
    }
    return Acts::Result<void>::failure(sResult.error());
  }
};

}  // namespace ActsFatras