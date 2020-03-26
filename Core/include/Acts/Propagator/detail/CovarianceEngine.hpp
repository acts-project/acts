// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class Surface;

namespace detail {

/// Construct bound parameters at the current position.
///
/// @param [in] freeParams free parameters vector
/// @param [in] boundCov bound covariance matrix
/// @param [in] covIsValid whether the covariance matrix contains valid entries
/// @param [in] surface surface on which the parameters are bound to
/// @param [in] geoCtx geometry context
/// @return bound parameters on the surface
///
/// @note Assumes that the position is already on the surface and covariance
///       (optionally) has already been transported.
BoundParameters makeBoundParameters(const FreeVector& freeParams,
                                    const BoundSymMatrix& boundCov,
                                    bool covIsValid, const Surface& surface,
                                    const GeometryContext& geoCtx);

/// Construct a curvilinear state at the current position
///
/// @param [in] freeParams free parameters vector
/// @param [in] curvilinearCov curvilinear covariance matrix
/// @param [in] covIsValid whether the covariance matrix contains valid entries
/// @return curvilinear parameters at the current position
///
/// @note Assumes that the covariance (optionally) has already been transported.
CurvilinearParameters makeCurvilinearParameters(
    const FreeVector& freeParams, const BoundSymMatrix& curvilinearCov,
    bool covIsValid);

/// Compute the combined jacobian to bound parameters on a surface.
///
/// @param [in] freeParams free parameters vector
/// @param [in] surface target surface for the bound parameters
/// @param [in] geoCtx geometry context
/// @param [in,out] jacToGlobal Jacobian from initial bound to free parameters
/// @param [in,out] jacTransport Transport Jacobian in free parameters
/// @param [in,out] derivative Stepper derivative
/// @param [out] jacobian from initial bound to target bound parameters.
///
/// @note This assumes that the current global parameter state has already
///       been propagated to be on the surface.
///
/// The partial Jacobians are reset such that they represent the propagation
/// starting at the surface.
void updateJacobiansToBound(const FreeVector& freeParams,
                            const Surface& surface,
                            const GeometryContext& geoCtx,
                            BoundToFreeMatrix& jacToGlobal,
                            FreeMatrix& jacTransport, FreeVector& derivative,
                            BoundMatrix& jacobian);

/// Compute the combined Jacobian to curvilinear parameters.
///
/// @param [in] freeParams free parameters vector
/// @param [in,out] jacToGlobal Jacobian from initial bound to free parameters
/// @param [in,out] jacTransport Transport Jacobian in free parameters
/// @param [in,out] derivative Stepper derivative
/// @param [out] jacobian from initial bound to target curvilinear parameters.
///
/// The partial Jacobians are reset such that they represent the propagation
/// starting at current curvilinear frame.
void updateJacobiansToCurvilinear(const FreeVector& freeParams,
                                  BoundToFreeMatrix& jacToGlobal,
                                  FreeMatrix& jacTransport,
                                  FreeVector& derivative,
                                  BoundMatrix& jacobian);

}  // namespace detail
}  // namespace Acts
