// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

class Surface;

namespace detail {

/// Convert global position and direction to bound parameters.
///
/// @param position Global track three-position
/// @param time Global track time
/// @param direction Global direction three-vector; normalization is ignored.
/// @param qOverP Charge-over-momentum-like parameter
/// @param surface Surface onto which the parameters are bound
/// @param geoCtx Geometry context for the global-to-local transformation
/// @return Bound parameters vector on the given surface
///
/// @warning The global position is assumed to be on the surface. If this is not
///          the case, the behaviour is undefined.
BoundVector transformFreeToBoundParameters(
    const Vector3D& position, double time, const Vector3D& direction,
    double qOverP, const Surface& surface, const GeometryContext& geoCtx);

/// Convert global direction to bound curvilinear parameters.
///
/// @param time Global track time
/// @param direction Global direction three-vector; normalization is ignored.
/// @param qOverP Charge-over-momentum-like parameter
/// @return Bound parameters vector on the curvilinear surface
///
/// @note The parameters are assumed to be defined at the origin of the
///       curvilinear frame derived from the direction vector. The local
///       coordinates are zero by construction.
BoundVector transformFreeToCurvilinearParameters(double time,
                                                 const Vector3D& direction,
                                                 double qOverP);

}  // namespace detail
}  // namespace Acts
