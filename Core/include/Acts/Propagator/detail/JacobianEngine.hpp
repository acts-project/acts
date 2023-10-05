// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// @brief These functions perform the calculation of the Jacobians for the
/// the covariance transport. This is a purely algebraic problem the
/// calculations are identical for @c StraightLineStepper and @c EigenStepper.
/// As a consequence the methods can be located in a separate file.
namespace detail {

/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
/// without transport jacobian.
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
FreeToBoundMatrix freeToCurvilinearJacobian(const Vector3& direction);

/// @brief Evaluate the projection Jacobian from curvilinear to free parameters
/// without transport jacobian.
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
BoundToFreeMatrix curvilinearToFreeJacobian(const Vector3& direction);

/// @brief This function calculates the jacobian from a free parameterisation
/// to a mixed one with angular representation.
///
/// @param direction The direction of the track
///
/// @return the 7x8 jacobian for direction to angles relation
ActsMatrix<7, 8> directionToAnglesJacobian(const Vector3& direction);

/// @brief This function calculates the jacobian from a free parameterisation
/// to a mixed one with angular representation.
///
/// @param direction The direction of the track
///
/// @return the 8x7 jacobian for angles to direction relation
ActsMatrix<8, 7> anglesToDirectionJacobian(const Vector3& direction);

/// @brief This function calculates the full transport jacobian from a bound
/// curvilinear representation to a new bound representation
///
/// @note Modifications of the jacobian related to the
/// projection onto a surface is considered. Since a variation of the start
/// parameters within a given uncertainty would lead to a variation of the end
/// parameters, these need to be propagated onto the target surface. This an
/// approximated approach to treat the (assumed) small change.
///
/// @param [in] geoContext The geometry Context
/// @param [in] freeParameters Free, nominal parametrisation
/// @param [in] boundToFreeJacobian Jacobian from bound to free at start
/// @param [in] freeTransportJacobian Transport jacobian free to free
/// @param [in] freeToPathDerivatives Path length derivatives for free
/// parameters
/// @param [in] surface Target surface
///
/// @note jac(locA->locB) = jac(gloB->locB)*(1+
/// pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
///
/// @return a 6x6 transport jacobian from bound to bound
BoundMatrix boundToBoundTransportJacobian(
    const GeometryContext& geoContext, const FreeVector& freeParameters,
    const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives, const Surface& surface);

/// @brief This function calculates the full jacobian from a given
/// bound/curvilinear parameterisation from a surface to new curvilinear
/// parameterisation.
///
/// @note Modifications of the jacobian related to the
/// projection onto a curvilinear surface is considered. Since a variation of
/// the start parameters within a given uncertainty would lead to a variation of
/// the end parameters, these need to be propagated onto the target surface.
/// This is an approximated approach to treat the (assumed) small change.
///
/// @param [in] direction Normalised direction vector
/// @param [in] boundToFreeJacobian Jacobian from bound to free at start
/// @param [in] freeTransportJacobian Transport jacobian free to free
/// @param [in] freeToPathDerivatives Path length derivatives for free
/// parameters
///
/// @note The parameter @p surface is only required if projected to bound
/// parameters. In the case of curvilinear parameters the geometry and the
/// position is known and the calculation can be simplified
///
/// @return a 6x6 transport jacobian from curilinear to bound
BoundMatrix boundToCurvilinearTransportJacobian(
    const Vector3& direction, const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives);

/// @brief This function calculates the full jacobian from a given
/// bound/curvilinear parameterisation from a new free
/// parameterisation.
///
/// @param [in] boundToFreeJacobian Jacobian from bound to free at start
/// @param [in] freeTransportJacobian Transport jacobian free to free
///
/// @return a 8x6 transport jacobian from bound to free
BoundToFreeMatrix boundToFreeTransportJacobian(
    const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian);

/// @brief This function calculates the full jacobian from a given
/// free parameterisation to a new curvilinear bound
/// parameterisation.
///
/// @note Modifications of the jacobian related to the
/// projection onto the target surface is considered. Since a variation of
/// the start parameters within a given uncertainty would lead to a variation of
/// the end parameters, these need to be propagated onto the target surface.
/// This is an approximated approach to treat the (assumed) small change.
///
/// @param [in] geoContext The geometry Context
/// @param [in] freeParameters Free, nominal parametrisation
/// @param [in] directionToAnglesJacobian The relation jacobian from dir to
/// angle
/// @param [in] anglesToDirectionJacobian The relation jacobian from angle to
/// dir
/// @param [in] freeTransportJacobian Transport jacobian free to free
/// @param [in] freeToPathDerivatives Path length derivatives for free
/// parameters
/// @param [in] surface The target surface
///
/// @return the 6x8 transport jacobian from bound to free
FreeToBoundMatrix freeToBoundTransportJacobian(
    const GeometryContext& geoContext, const FreeVector& freeParameters,
    const ActsMatrix<7, 8>& directionToAnglesJacobian,
    const ActsMatrix<8, 7>& anglesToDirectionJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives, const Surface& surface);

/// @brief This function calculates the full transport jacobian from a free
/// parameterisation to a bound one. Since a variation of the start
/// parameters within a given uncertainty would lead to a variation of the end
/// parameters, these need to be propagated onto the target surface. This an
/// approximated approach to treat the (assumed) small change.
///
/// @param [in] direction Normalised direction vector
/// @param [in] directionToAnglesJacobian The relation jacobian from dir to
/// angle
/// @param [in] anglesToDirectionJacobian The relation jacobian from angle to
/// dir
/// @param [in] freeTransportJacobian Transport jacobian free to free
/// @param [in] freeToPathDerivatives Path length derivatives for free
///
/// @return a 6x8 transport jacobian from curvilinear to free
FreeToBoundMatrix freeToCurvilinearTransportJacobian(
    const Vector3& direction, const ActsMatrix<7, 8>& directionToAnglesJacobian,
    const ActsMatrix<8, 7>& anglesToDirectionJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives);

/// @brief This function calculates the free transfport jacobian from a free
/// parameterisation.
/// @param [in] directionToAnglesJacobian The relation jacobian from dir to
/// angle
/// @param [in] anglesToDirectionJacobian The relation jacobian from angle to
/// dir
/// @param [in] freeTransportJacobian Transport jacobian free to free
///
/// @return a 8x8 transport jacobian from free to free
FreeMatrix freeToFreeTransportJacobian(
    const ActsMatrix<7, 8>& directionToAnglesJacobian,
    const ActsMatrix<8, 7>& anglesToDirectionJacobian,
    const FreeMatrix& freeTransportJacobian);

}  // namespace detail

}  // namespace Acts
