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

/// @brief Evaluate the derivative of bound parameters w.r.t. alignment
/// parameters (i.e. position and rotation) of its reference surface
///
/// @param [in] geoContext The geometry Context
/// @param [in] boundParams The bound parameters to investigate
/// @param [in] derivatives Path length derivatives of the free, nominal
/// parameters to help evaluate path correction due to change of alignment
/// parameters
///
/// @return Derivative of bound parameters w.r.t. alignment parameters
AlignmentToBoundMatrix alignmentToLocalDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives) const;

}  // namespace detail
}  // namespace Acts
