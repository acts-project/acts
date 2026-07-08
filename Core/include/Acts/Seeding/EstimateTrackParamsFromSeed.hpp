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
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <system_error>
#include <type_traits>

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

/// Estimate free track parameters from three space points.
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

/// Estimate free track parameters from three space points.
///
/// Optionally estimates the tangents at the three space points from helix fit.
/// This can be used as an input for the strip space point calibration, which
/// requires the track tangents at each space point as input.
///
/// @param sp0 is the bottom space point
/// @param t0 is the time of the bottom space point
/// @param sp1 is the middle space point
/// @param sp2 is the top space point
/// @param bField is the magnetic field vector
/// @param tangent0 is the output tangent at the bottom space point
/// @param tangent1 is the output tangent at the middle space point
/// @param tangent2 is the output tangent at the top space point
///
/// @return the free parameters
FreeVector estimateTrackParamsFromSeed(const Vector3& sp0, double t0,
                                       const Vector3& sp1, const Vector3& sp2,
                                       const Vector3& bField,
                                       Vector3* tangent0 = nullptr,
                                       Vector3* tangent1 = nullptr,
                                       Vector3* tangent2 = nullptr);

/// Estimate free track parameters from three space points
///
/// @tparam space_point_range_t The type of space point range
///
/// @param spRange is the range of space points
/// @param bField is the magnetic field vector
///
/// @return the free parameters
template <std::ranges::range space_point_range_t>
[[deprecated(
    "The broadly templated versions of estimateTrackParamsFromSeed will be "
    "removed in the future.")]]
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
[[deprecated(
    "The broadly templated versions of estimateTrackParamsFromSeed will be "
    "removed in the future.")]]
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

/// Error codes for the multi-space-point track parameter estimation
/// @ingroup errors
enum class TrackParamsEstimationError {
  // ensure all values are non-zero
  /// Fewer than three space points were provided
  NotEnoughSpacePoints = 1,
  /// The fit is degenerate (e.g. all space points coincide)
  DegenerateFit,
};

/// Create error code from @ref TrackParamsEstimationError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(Acts::TrackParamsEstimationError e);

/// Estimate free track parameters from an ordered set of N >= 3 space points.
///
/// Generalizes @ref estimateTrackParamsFromSeed (limited to three space points)
/// to a least-squares helix fit: the circle transverse to the field is fitted
/// algebraically with the Taubin method (optionally refined geometrically) and
/// the coordinate along the field is fitted linearly against the transverse arc
/// length. The parameters are represented at the first (reference) space point.
/// Points are taken in track order (first = reference, usually innermost) and
/// are not sorted internally.
///
/// Straight tracks are handled naturally: for a vanishing curvature (high
/// momentum or weak field) the circle fit degenerates to a line and only the
/// direction is estimated. q/p is left at zero when the field magnitude
/// vanishes, as the momentum cannot be estimated without a field.
///
/// An optional per-space-point weight turns every least-squares stage (the
/// circle fit, the R-Z line fit and the straight-line PCA fallback) into a
/// weighted fit, expressing the relative trust in each point (e.g. a more
/// precise detector type). A single scalar weight per point is used rather than
/// a per-plane split: the estimation frame is a general rotation fixed by the
/// field, so a transverse/longitudinal variance split would not be
/// rotation-invariant and is only meaningful when the field is along global z.
/// Weights act as relative (e.g. inverse-variance) factors; an empty span
/// selects uniform weights, so passing none reproduces the unweighted fit. A
/// non-empty span must match `spacePoints` in size.
///
/// @param spacePoints the ordered global space point positions
/// @param bField the homogeneous magnetic field vector
/// @param t0 the time assigned to the reference point (eFreeTime)
/// @param geometricRefineIterations number of Gauss-Newton refinement
///        iterations on top of the algebraic circle fit (0 disables it)
/// @param weights optional per-point weights for all fit stages
///        (empty span = uniform)
/// @return the free parameters at the reference point, or an error
Result<FreeVector> estimateTrackParamsFromSpacePoints(
    std::span<const Vector3> spacePoints, const Vector3& bField, double t0 = 0.,
    std::size_t geometricRefineIterations = 0,
    std::span<const double> weights = {});

/// Estimate bound track parameters from an ordered set of N >= 3 space points.
///
/// As @ref estimateTrackParamsFromSpacePoints, but expressed at the given
/// surface, which is assumed to be that of the first (reference) space point.
///
/// @param gctx the geometry context
/// @param surface the surface of the reference space point
/// @param spacePoints the ordered global space point positions
/// @param bField the homogeneous magnetic field vector
/// @param t0 the time assigned to the reference point (eBoundTime)
/// @param geometricRefineIterations number of Gauss-Newton refinement
///        iterations on top of the algebraic circle fit (0 disables it)
/// @param weights optional per-point weights for all fit stages
///        (empty span = uniform)
/// @return the bound parameters at the surface, or an error
Result<BoundVector> estimateTrackParamsFromSpacePoints(
    const GeometryContext& gctx, const Surface& surface,
    std::span<const Vector3> spacePoints, const Vector3& bField, double t0 = 0.,
    std::size_t geometricRefineIterations = 0,
    std::span<const double> weights = {});

/// @}

}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::TrackParamsEstimationError> : std::true_type {};
}  // namespace std
