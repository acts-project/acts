// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
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

/// Transform bound track parameters into equivalent free track parameters.
///
/// @param surface Surface onto which the input parameters are bound
/// @param geoCtx Geometry context for the local-to-global transformation
/// @param boundParams Bound track parameters vector
/// @return Equivalent free trackparameters vector
FreeVector transformBoundToFreeParameters(const Surface& surface,
                                          const GeometryContext& geoCtx,
                                          const BoundVector& boundParams);

}  // namespace detail
}  // namespace Acts
