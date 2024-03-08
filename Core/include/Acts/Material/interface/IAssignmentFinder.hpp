// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialInteraction.hpp"

#include <tuple>
#include <utility>
#include <vector>

#pragma once

namespace Acts {

/// @brief Interface for the material mapping that seeks the possibile
/// assignment candidates for the material interactiosn
class IAssignmentFinder {
 public:
  /// Virtual destructor
  virtual ~IAssignmentFinder() = default;

  /// @brief SurfaceAssignment is a surface, a position and a direction
  using SurfaceAssignment = std::tuple<const Surface*, Vector3, Vector3>;

  /// @brief VolumeAssignment is a volume and a entry and exit of the volume
  using VolumeAssignment =
      std::tuple<const InteractionVolume, Vector3, Vector3>;

  /// @brief Interface method for generating assigment candidates for the
  /// material interaction assignment to surfaces or volumes
  ///
  /// @param gctx is the geometry context
  /// @param mctx is the magnetic field context
  /// @param position is the position of the initial ray
  /// @param direction is the direction of initial ray
  ///
  /// @return a vector of Surface Assignments and Volume Assignments
  virtual std::pair<std::vector<SurfaceAssignment>,
                    std::vector<VolumeAssignment>>
  assignmentCandidates(const GeometryContext& gctx,
                       const MagneticFieldContext& mctx,
                       const Vector3& position,
                       const Vector3& direction) const = 0;
};

}  // namespace Acts
