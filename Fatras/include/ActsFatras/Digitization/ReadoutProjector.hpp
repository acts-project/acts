// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Digitization/DigitizationData.hpp"
#include "ActsFatras/Digitization/detail/PlanarSurfaceMask.hpp"
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Result.hpp>

namespace ActsFatras {

using ProjectedPosition = std::pair<Acts::Vector2D, double>;
using ProjectedSegment = std::array<ProjectedPosition, 2>;

/// Project an impact point onto the surface including drift direction
///
struct ReadoutProjector {
  /// The projection onto the readout surface.
  ///
  /// @param dInput The digitization input including the hit and surface
  /// @param dDir The drift direction (check with w.r.t
  /// @param dThickness The module thickness
  ///
  /// @return The drifted position
  Acts::Result<ProjectedSegment> project(const DigitizationInput& dInput,
                                         const Acts::Vector3D& dDir,
                                         double dThickness) const;

  /// The surface maksdre
  detail::PlanarSurfaceMask surfaceMask;
};

}  // namespace ActsFatras