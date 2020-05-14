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
using AlignmentToCartesian3DMatrix =
    ActsMatrix<BoundParametersScalar, 3, eAlignmentParametersSize>;
using Cartesian3DToLocal2DMatrix = ActsMatrix<BoundParametersScalar, 2, 3>;

/// @brief Evaluate the derivative of bound track parameters w.r.t. alignment
/// parameters (i.e. local reference frame origin in global 3D cartesian
/// coordinates and rotation represented with Euler angles)
///
/// @param [in] geoContext The geometry Context
/// @param [in] boundParams The bound parameters to investigate
/// @param [in] derivatives Path length derivatives of the free, nominal
/// parameters to help evaluate path correction due to change of alignment
/// parameters
/// @param [in] rframeOrigin The origin of local reference frame in global
/// coordinate
/// @param [in] cartesianToLocal The derivative of track position represented in
/// (local) bound track parameters (could be in non-cartesian coordinates)
/// w.r.t. track position represented in local 3D cartesian coordinates. This is
/// needed because alignment is done w.r.t. cartesian coordinates
///
/// @return Derivative of bound track parameters w.r.t. alignment parameters
AlignmentToBoundMatrix alignmentToLocalDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives, const Vector3D& rframeOrigin,
    const Cartesian3DToLocal2DMatrix& cartesianToLocal) const;

}  // namespace detail
}  // namespace Acts
