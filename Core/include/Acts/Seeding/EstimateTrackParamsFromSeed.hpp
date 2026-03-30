// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <array>
#include <optional>
#include <stdexcept>

namespace Acts {

/// @defgroup est_track_params Estimate track parameters from seed
///
/// The implemented method is based on the conformal map transformation. It
/// estimates the full free track parameters, i.e. (x, y, z, t, dx, dy, dz, q/p)
/// at the bottom space point. The magnetic field (which can be along any
/// direction) is also necessary for the momentum estimation.
///
/// It resembles the method used in ATLAS for the track parameters estimated
/// from seed, i.e. the function InDet::SiTrackMaker_xk::getAtaPlane here:
/// https://acode-browser.usatlas.bnl.gov/lxr/source/athena/InnerGeometry/InDetRecTools/SiTrackMakerTool_xk/src/SiTrackMaker_xk.cxx
///
/// @{

/// Estimate free track parameters from three space points
///
/// This is a purely spatial estimation, i.e. the time parameter will be set to
/// 0.
///
/// @param sp0 is the bottom space point
/// @param sp1 is the middle space point
/// @param sp2 is the top space point
/// @param bField is the magnetic field vector
///
/// @return the free parameters
[[deprecated(
    "Use the version of estimateTrackParamsFromSeed with time information "
    "instead.")]]
FreeVector estimateTrackParamsFromSeed(const Vector3& sp0, const Vector3& sp1,
                                       const Vector3& sp2,
                                       const Vector3& bField);

/// Estimate free track parameters from three space points
///
/// @param sp0 is the bottom space point
/// @param t0 is the time of the bottom space point
/// @param sp1 is the middle space point
/// @param sp2 is the top space point
/// @param bField is the magnetic field vector
///
/// @return the free parameters
FreeVector estimateTrackParamsFromSeed(const Vector3& sp0, double t0,
                                       const Vector3& sp1, const Vector3& sp2,
                                       const Vector3& bField);

/// Estimate free track parameters from three space points
///
/// @tparam space_point_range_t The type of space point range
///
/// @param spRange is the range of space points
/// @param bField is the magnetic field vector
///
/// @return the free parameters
template <std::ranges::range space_point_range_t>
FreeVector estimateTrackParamsFromSeed(space_point_range_t spRange,
                                       const Vector3& bField) {
  // Check the number of provided space points
  if (spRange.size() != 3) {
    throw std::invalid_argument(
        "There should be exactly three space points provided.");
  }

  // The global positions of the bottom, middle and space points
  std::array<Vector3, 3> spPositions = {Vector3::Zero(), Vector3::Zero(),
                                        Vector3::Zero()};
  std::array<std::optional<double>, 3> spTimes = {std::nullopt, std::nullopt,
                                                  std::nullopt};
  // The first, second and third space point are assumed to be bottom, middle
  // and top space point, respectively
  for (auto [sp, spPosition, spTime] :
       Acts::zip(spRange, spPositions, spTimes)) {
    if (sp == nullptr) {
      throw std::invalid_argument("Empty space point found.");
    }
    spPosition = Vector3(sp->x(), sp->y(), sp->z());
    spTime = sp->t();
  }

  return estimateTrackParamsFromSeed(spPositions[0], spTimes[0].value_or(0),
                                     spPositions[1], spPositions[2], bField);
}

/// Estimate bound track parameters from three space points
///
/// @param gctx is the geometry context
/// @param surface is the surface of the bottom space point. The estimated bound
///                track parameters will be represented at this surface.
/// @param sp0 is the bottom space point
/// @param t0 is the time of the bottom space point
/// @param sp1 is the middle space point
/// @param sp2 is the top space point
/// @param bField is the magnetic field vector
///
/// @return bound parameters
Result<BoundVector> estimateTrackParamsFromSeed(
    const GeometryContext& gctx, const Surface& surface, const Vector3& sp0,
    double t0, const Vector3& sp1, const Vector3& sp2, const Vector3& bField);

/// Estimate bound track parameters from three space points
///
/// @tparam space_point_range_t The type of space point range
///
/// @param gctx is the geometry context
/// @param spRange is the range of space points
/// @param surface is the surface of the bottom space point. The estimated bound
///                track parameters will be represented at this surface.
/// @param bField is the magnetic field vector
///
/// @return bound parameters
template <std::ranges::range space_point_range_t>
Result<BoundVector> estimateTrackParamsFromSeed(const GeometryContext& gctx,
                                                space_point_range_t spRange,
                                                const Surface& surface,
                                                const Vector3& bField) {
  const FreeVector freeParams = estimateTrackParamsFromSeed(spRange, bField);
  return transformFreeToBoundParameters(freeParams, surface, gctx);
}

/// Configuration for the estimation of the covariance matrix of the track
/// parameters with `estimateTrackParamCovariance`.
struct EstimateTrackParamCovarianceConfig {
  /// The initial sigmas for the track parameters
  BoundVector initialSigmas = {1. * UnitConstants::mm,
                               1. * UnitConstants::mm,
                               1. * UnitConstants::degree,
                               1. * UnitConstants::degree,
                               1. * UnitConstants::e / UnitConstants::GeV,
                               1. * UnitConstants::ns};

  /// The initial sigma for the q/pt
  /// @note The resulting q/p sigma is added to the one in `initialSigmas`
  double initialSigmaQoverPt = 0. * UnitConstants::e / UnitConstants::GeV;

  /// The initial relative uncertainty sigma(pt)/pt
  /// @note The resulting q/p sigma is added to the one in `initialSigmas`
  double initialSigmaPtRel = 0.1;

  /// The inflation factors for the variances of the track parameters
  BoundVector initialVarInflation = {1., 1., 1., 1., 1., 1.};
  /// The inflation factor for time uncertainty if the time parameter was not
  /// estimated
  double noTimeVarInflation = 100.;
};

/// Estimate the covariance matrix of the given track parameters based on the
/// provided configuration. The assumption is that we can model the uncertainty
/// of the track parameters as a diagonal matrix with the provided initial
/// sigmas. The inflation factors are used to inflate the initial variances
/// based on the provided configuration. The uncertainty of q/p is estimated
/// based on the relative uncertainty of the q/pt and the theta uncertainty.
///
/// @param config is the configuration for the estimation
/// @param params is the track parameters
/// @param hasTime is true if the track parameters have time
///
/// @return the covariance matrix of the track parameters
BoundMatrix estimateTrackParamCovariance(
    const EstimateTrackParamCovarianceConfig& config, const BoundVector& params,
    bool hasTime);

/// @}

}  // namespace Acts
