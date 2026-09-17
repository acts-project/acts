// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"

#include <optional>

namespace Acts {
class Surface;
}

namespace Acts::detail::sympy {

/// @brief These functions perform the transport of a covariance matrix using
/// given Jacobians. The required data is provided by the stepper object
/// with some additional data. Since this is a purely algebraic problem the
/// calculations are identical for @c StraightLineStepper and @c EigenStepper.
/// As a consequence the methods can be located in a separate file.

/// @brief Method for on-demand covariance transport of a bound/curvilinear to
///        another bound representation.
///
/// @param [in] geoContext The geometry context
/// @param [in] surface is the surface to which the covariance is forwarded to
/// @param [in, out] boundCovariance The covariance matrix of the state
/// @param [in, out] fullTransportJacobian Full jacobian since the last reset
/// @param [in, out] freeToPathDerivatives Path length derivatives
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
///        parametrisation to free parameters
/// @param [in, out] freeParameters Free, nominal parametrisation
/// @param [in] freeToBoundCorrection Correction for non-linearity effect during
///        transform from free to bound
///
/// @note No check is done if the position is actually on the surface
///
/// @return Failure if the parameters cannot be expressed on the surface
Result<void> transportCovarianceToBound(
    const GeometryContext& geoContext, const Surface& surface,
    BoundMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeVector& freeToPathDerivatives, BoundToFreeMatrix& boundToFreeJacobian,
    const std::optional<FreeMatrix>& additionalFreeCovariance,
    FreeVector& freeParameters,
    const FreeToBoundCorrection& freeToBoundCorrection);

/// @brief Method for on-demand covariance transport of a bound/curvilinear
///        to a new curvilinear representation.
///
/// @param [in, out] boundCovariance The covariance matrix of the state
/// @param [in, out] fullTransportJacobian Full jacobian since the last reset
/// @param [in, out] freeToPathDerivatives Path length derivatives
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
///        parametrisation to free parameters
/// @param [in] direction Normalised direction vector
///
void transportCovarianceToCurvilinear(
    BoundMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeVector& freeToPathDerivatives, BoundToFreeMatrix& boundToFreeJacobian,
    const std::optional<FreeMatrix>& additionalFreeCovariance,
    const Vector3& direction);

}  // namespace Acts::detail::sympy
