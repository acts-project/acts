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
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"

#include <optional>

namespace Acts::detail {

/// @brief These functions perform the transport of a covariance matrix using
/// given Jacobians. The required data is provided by the stepper object
/// with some additional data. Since this is a purely algebraic problem the
/// calculations are identical for @c StraightLineStepper and @c EigenStepper.
/// As a consequence the methods can be located in a separate file.

/// Create bound parameters from free parameters on a surface
///
/// @note It does not check if the free parameters are on the surface
///
/// @param [in] geoContext The geometry context
/// @param [in] surface The surface of the bound parameters
/// @param [in] freeParameters The free parameters
/// @param [in] covariance The optional bound covariance on @p surface
/// @param [in] particleHypothesis The particle hypothesis
///
/// @return The bound parameters, or a failure if the free parameters cannot
///         be expressed on @p surface
Result<BoundTrackParameters> boundParameters(
    const GeometryContext& geoContext, const Surface& surface,
    const FreeVector& freeParameters, std::optional<BoundMatrix> covariance,
    const ParticleHypothesis& particleHypothesis);

/// Create curvilinear parameters from free parameters
///
/// @param [in] freeParameters The free parameters
/// @param [in] covariance The optional curvilinear covariance
/// @param [in] particleHypothesis The particle hypothesis
///
/// @return The curvilinear parameters at the free position
BoundTrackParameters curvilinearParameters(
    const FreeVector& freeParameters, std::optional<BoundMatrix> covariance,
    const ParticleHypothesis& particleHypothesis);

/// @brief Method for on-demand covariance transport of a bound/curvilinear to
///        another bound representation.
///
/// @param [in] geoContext The geometry context
/// @param [in] surface is the surface to which the covariance is forwarded to
/// @param [in, out] boundCovariance The covariance matrix of the state
/// @param [in, out] fullTransportJacobian Full jacobian since the last reset
/// @param [in, out] freeTransportJacobian Global jacobian since the last reset
/// @param [in, out] freeToPathDerivatives Path length derivatives
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
///        parametrisation to free parameters
/// @param [in, out] freeParameters Free, nominal parametrisation
/// @param [in] freeToBoundCorrection Correction for non-linearity effect during
///        transform from free to bound
///
/// @return Failure if the parameters are not on the surface
Result<void> transportCovarianceToBound(
    const GeometryContext& geoContext, const Surface& surface,
    BoundMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian,
    const std::optional<FreeMatrix>& additionalFreeCovariance,
    FreeVector& freeParameters,
    const FreeToBoundCorrection& freeToBoundCorrection);

/// @brief Method for on-demand covariance transport of a bound/curvilinear
///        to a new curvilinear representation.
///
/// @param [in, out] boundCovariance The covariance matrix of the state
/// @param [in, out] fullTransportJacobian Full jacobian since the last reset
/// @param [in, out] freeTransportJacobian Global jacobian since the last reset
/// @param [in, out] freeToPathDerivatives Path length derivatives
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
///        parametrisation to free parameters
/// @param [in] direction Normalised direction vector
///
void transportCovarianceToCurvilinear(
    BoundMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian,
    const std::optional<FreeMatrix>& additionalFreeCovariance,
    const Vector3& direction);

/// Convert bound track parameters to another bound surface.
/// @pre The @p targetSurface must intersect with the surface attached to
///      @p boundParameters, and the parameters must be on-surface on the
///      target surface.
/// @param gctx The geometry context.
/// @param boundParameters The bound track parameters to convert.
/// @param targetSurface The target surface.
/// @param bField The magnetic field at the target surface.
/// @return The converted bound track parameters.
Result<BoundTrackParameters> boundToBoundConversion(
    const GeometryContext& gctx, const BoundTrackParameters& boundParameters,
    const Surface& targetSurface, const Vector3& bField = Vector3::Zero());

}  // namespace Acts::detail
