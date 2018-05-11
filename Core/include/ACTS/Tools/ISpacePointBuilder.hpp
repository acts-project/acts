// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_TOOLS_ISPACEPOINTBUILDER_H
#define ACTS_TOOLS_ISPACEPOINTBUILDER_H

#include <iterator>
#include <vector>
#include "ACTS/Digitization/PlanarModuleCluster.hpp"

namespace Acts {

/// @brief Structure for easier bookkeeping of hits.
struct SpacePoint
{
  /// Storage of the hit(s) on a surface/surfaces
  std::vector<PlanarModuleCluster const*> hitModule;
  /// Storage of a space point. Zero vector indicates unset point
  Vector3D spacePoint = {0., 0., 0.};
};

/// @class ISpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the hits on some detector elements need further treatment. This interface
/// serves to take the digitized hits on a detector element and provide the
/// corresponding space point.
///
class ISpacePointBuilder
{
public:
  /// @brief Adds hits to the list. Allows checks for possible combination
  /// vetos.
  /// @param spacePoint storage of the space points
  /// @param hits list of list of hits. The 2D setup allows possible combination
  /// vetos
  virtual void
  addHits(std::vector<SpacePoint>& spacePoints,
          const std::vector<std::vector<Acts::PlanarModuleCluster const*>>&
              hits) const = 0;

  /// @brief Calculates the space points out of a given collection of hits and
  /// stores the data
  /// @param spacePoints storage of the data
  virtual void
  calculateSpacePoints(std::vector<SpacePoint>& spacePoints) const = 0;

protected:
  /// @brief Getter method for the local coordinates of a hit
  /// on its corresponding surface
  /// @param hit object related to the hit that holds the necessary information
  /// @return vector of the local coordinates of the hit on the surface
  virtual Vector2D
  localCoords(const PlanarModuleCluster& hit) const = 0;

  /// @brief Getter method for the global coordinates of a hit
  /// @param hit object related to the hit that holds the necessary information
  /// @return vector of the global coordinates of the hit
  virtual Vector3D
  globalCoords(const PlanarModuleCluster& hit) const = 0;
};

}  // namespace Acts

#endif  // ACTS_TOOLS_ISPACEPOINTFINDER_H
