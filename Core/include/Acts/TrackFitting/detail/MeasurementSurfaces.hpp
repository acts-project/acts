// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts::detail {

/// Make the navigator reach a measurement surface without a bounds check. A
/// surface of the tracking geometry is extended. A free surface has no
/// geometry identifier, and the navigator throws for it as an extended
/// surface, so it is offered as an additional surface.
/// @param options The navigator options of the fit
/// @param surface The measurement surface
inline void registerMeasurementSurface(NavigatorPlainOptions& options,
                                       const Surface& surface) {
  if (surface.geometryId() == GeometryIdentifier{}) {
    options.registerAdditionalSurface(surface);
  } else {
    options.registerExtendedSurface(surface);
  }
}

}  // namespace Acts::detail
