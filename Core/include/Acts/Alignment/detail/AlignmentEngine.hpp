// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace detail {
using AlignmentToCartesianMatrix =
    ActsMatrix<AlignmentParametersScalar, eCartesianCoordinatesDimension,
               eAlignmentParametersSize>;
using CartesianToBoundLocalMatrix =
    ActsMatrix<BoundParametersScalar, 2, eCartesianCoordinatesDimension>;

/// @brief Evaluate the derivative of bound track parameters w.r.t. alignment
/// parameters of its reference surface (i.e. local reference frame origin in
/// global 3D Cartesian coordinates and rotation represented with extrinsic
/// Euler angles)
///
/// @param [in] geoContext The geometry Context
/// @param [in] boundParams The bound parameters to investigate
/// @param [in] derivatives Path length derivatives of the free, nominal
/// parameters to help evaluate path correction due to change of alignment
/// parameters
/// @param [in] rframeOrigin The origin of local reference frame in global
/// coordinate
/// @param [in] cartesianToLocal The derivative of track position represented in
/// (local) bound track parameters (could be in non-Cartesian coordinates)
/// w.r.t. track position represented in local 3D Cartesian coordinates. This is
/// needed because alignment is done w.r.t. Cartesian coordinates
///
/// @return Derivative of bound track parameters w.r.t. local reference frame
/// alignment parameters
AlignmentToBoundMatrix surfaceAlignmentToBoundDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives, const Vector3D& rframeOrigin,
    const CartesianToBoundLocalMatrix& locCartesianToLocBound =
        CartesianToBoundLocalMatrix::Identity(2,
                                              eCartesianCoordinatesDimension));

/// @brief Evaluate the derivative of bound track parameters w.r.t. alignment
/// parameters of its referece surface's associated layer (i.e. layer center in
/// global 3D Cartesian coordinates and layer rotation represented with
/// extrinsic Euler angles)
/// @Todo: The jacobian of layer alignment parameters to surface alignment
/// parameters could be calculated inside this function?
///
/// @param [in] geoContext The geometry Context
/// @param [in] boundParams The bound parameters to investigate
/// @param [in] derivatives Path length derivatives of the free, nominal
/// parameters to help evaluate path correction due to change of alignment
/// parameters
/// @param [in] rframeOrigin The origin of local reference frame in global
/// coordinate
/// @param [in] cartesianToLocal The derivative of track position represented in
/// (local) bound track parameters (could be in non-Cartesian coordinates)
/// w.r.t. track position represented in local 3D Cartesian coordinates. This is
/// needed because alignment is done w.r.t. Cartesian coordinates
/// @param [in] layerAlignToSurfaceAlign The jacobian of layer alignment
/// parameters to surface alignment parameters
///
/// @return Derivative of bound track parameters w.r.t. layer alignment
/// parameters
AlignmentToBoundMatrix layerAlignmentToBoundDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives, const Vector3D& rframeOrigin,
    const AlignmentMatrix& layerAlignToSurfaceAlign,
    const CartesianToBoundLocalMatrix& locCartesianToLocBound =
        CartesianToBoundLocalMatrix::Identity(2,
                                              eCartesianCoordinatesDimension));

/// @brief Evaluate the derivative of bound track parameters w.r.t. alignment
/// parameters of its referece surface's associated volume (i.e. volume center
/// in global 3D Cartesian coordinates and volume rotation represented with
/// extrinsic Euler angles)
/// @Todo: The jacobian of volume alignment parameters to surface alignment
/// parameters could be calculated inside this function?
///
/// @param [in] geoContext The geometry Context
/// @param [in] boundParams The bound parameters to investigate
/// @param [in] derivatives Path length derivatives of the free, nominal
/// parameters to help evaluate path correction due to change of alignment
/// parameters
/// @param [in] rframeOrigin The origin of local reference frame in global
/// coordinate
/// @param [in] cartesianToLocal The derivative of track position represented in
/// (local) bound track parameters (could be in non-Cartesian coordinates)
/// w.r.t. track position represented in local 3D Cartesian coordinates. This is
/// needed because alignment is done w.r.t. Cartesian coordinates
/// @param [in] volumeAlignToSurfaceAlign The jacobian of layer alignment
/// parameters to volume alignment parameters
///
/// @return Derivative of bound track parameters w.r.t. volume alignment
/// parameters
AlignmentToBoundMatrix volumeAlignmentToBoundDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives, const Vector3D& rframeOrigin,
    const AlignmentMatrix& volumeAlignToSurfaceAlign,
    const CartesianToBoundLocalMatrix& locCartesianToLocBound =
        CartesianToBoundLocalMatrix::Identity(2,
                                              eCartesianCoordinatesDimension));

/// @brief Evaluate the derivative of local reference frame axes vector w.r.t.
/// its rotation around global x/y/z axis
/// @Todo: add parameter for rotation axis order
///
/// @param [in] rframe The local reference frame
///
/// @return Derivative of local reference frame x/y/z axis vector w.r.t. its
/// rotation angles (extrinsic Euler angles) around global x/y/z axis
std::tuple<RotationMatrix3D, RotationMatrix3D, RotationMatrix3D>
rotationToLocalAxesDerivative(const RotationMatrix3D& rframe);

}  // namespace detail
}  // namespace Acts
