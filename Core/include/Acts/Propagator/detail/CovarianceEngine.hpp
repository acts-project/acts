// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <functional>
#include <tuple>
#include <variant>

namespace Acts {

/// @brief These functions perform the transport of a covariance matrix using
/// given Jacobians. The required data is provided by the stepper object
/// with some additional data. Since this is a purely algebraic problem the
/// calculations are identical for @c StraightLineStepper and @c EigenStepper.
/// As a consequence the methods can be located in a separate file.
namespace detail {

/// Create and return the bound state at the current position
///
/// @brief It does not check if the transported state is at the surface, this
/// needs to be guaranteed by the propagator
///
/// @param [in] geoContext The geometry context
/// @param [in, out] covarianceMatrix The covariance matrix of the state
/// @param [in, out] jacobian Full jacobian since the last reset
/// @param [in, out] transportJacobian Global jacobian since the last reset
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] jacToGlobal Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in, out] parameters Free, nominal parametrisation
/// @param [in] particleHypothesis Particle hypothesis
/// @param [in] covTransport Decision whether the covariance transport should be
/// performed
/// @param [in] accumulatedPath Propagated distance
/// @param [in] surface Target surface on which the state is represented
/// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
///
/// @return A bound state:
///   - the parameters at the surface
///   - the stepwise jacobian towards it (from last bound)
///   - and the path length (from start - for ordering)
Result<std::tuple<BoundTrackParameters, BoundMatrix, double>> boundState(
    const GeometryContext& geoContext, BoundSquareMatrix& covarianceMatrix,
    BoundMatrix& jacobian, FreeMatrix& transportJacobian,
    FreeVector& derivatives, BoundToFreeMatrix& jacToGlobal,
    FreeVector& parameters, const ParticleHypothesis& particleHypothesis,
    bool covTransport, double accumulatedPath, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection =
        FreeToBoundCorrection(false));

/// Create and return a curvilinear state at the current position
///
/// @brief This creates a curvilinear state.
///
/// @param [in, out] covarianceMatrix The covariance matrix of the state
/// @param [in, out] jacobian Full jacobian since the last reset
/// @param [in, out] transportJacobian Global jacobian since the last reset
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] jacToGlobal Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] parameters Free, nominal parametrisation
/// @param [in] particleHypothesis Particle hypothesis
/// @param [in] covTransport Decision whether the covariance transport should be
/// performed
/// @param [in] accumulatedPath Propagated distance
///
/// @return A curvilinear state:
///   - the curvilinear parameters at given position
///   - the stepweise jacobian towards it (from last bound)
///   - and the path length (from start - for ordering)
std::tuple<CurvilinearTrackParameters, BoundMatrix, double> curvilinearState(
    BoundSquareMatrix& covarianceMatrix, BoundMatrix& jacobian,
    FreeMatrix& transportJacobian, FreeVector& derivatives,
    BoundToFreeMatrix& jacToGlobal, const FreeVector& parameters,
    const ParticleHypothesis& particleHypothesis, bool covTransport,
    double accumulatedPath);

/// @brief Method for on-demand covariance transport of a bound/curvilinear to
/// another bound representation.
///
/// @param [in] geoContext The geometry context
/// @param [in, out] boundCovariance The covariance matrix of the state
/// @param [in, out] fullTransportJacobian Full jacobian since the last reset
/// @param [in, out] freeTransportJacobian Global jacobian since the last reset
/// @param [in, out] freeToPathDerivatives Path length derivatives
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in, out] freeParameters Free, nominal parametrisation
/// @param [in] surface is the surface to which the covariance is
///        forwarded to
/// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
///
/// @note No check is done if the position is actually on the surface
///
void transportCovarianceToBound(
    const GeometryContext& geoContext, BoundSquareMatrix& boundCovariance,
    BoundMatrix& fullTransportJacobian, FreeMatrix& freeTransportJacobian,
    FreeVector& freeToPathDerivatives, BoundToFreeMatrix& boundToFreeJacobian,
    FreeVector& freeParameters, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection =
        FreeToBoundCorrection(false));

/// @brief Method for on-demand covariance transport of a bound/curvilinear
/// to a new curvilinear representation.
///
/// @param [in, out] boundCovariance The covariance matrix of the state
/// @param [in, out] fullTransportJacobian Full jacobian since the last reset
/// @param [in, out] freeTransportJacobian Global jacobian since the last reset
/// @param [in, out] freeToPathDerivatives Path length derivatives
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] direction Normalised direction vector
///
void transportCovarianceToCurvilinear(BoundSquareMatrix& boundCovariance,
                                      BoundMatrix& fullTransportJacobian,
                                      FreeMatrix& freeTransportJacobian,
                                      FreeVector& freeToPathDerivatives,
                                      BoundToFreeMatrix& boundToFreeJacobian,
                                      const Vector3& direction);

}  // namespace detail
}  // namespace Acts
