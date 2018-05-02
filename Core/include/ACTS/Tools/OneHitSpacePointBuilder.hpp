// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_TOOLS_ONEHITSPACEPOINTBUILDER_H
#define ACTS_TOOLS_ONEHITSPACEPOINTBUILDER_H

#include "ACTS/Digitization/CartesianSegmentation.hpp"
#include "ACTS/Digitization/DigitizationModule.hpp"
#include "ACTS/Tools/ISpacePointBuilder.hpp"

namespace Acts {

/// @class OneHitSpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the hits strip detectors need further treatment. This class takes
/// the digitized hits on a strip or pixel detector element and provides the
/// corresponding space point.
///
class OneHitSpacePointBuilder : public ISpacePointBuilder
{
public:
  /// Constructor
  OneHitSpacePointBuilder();

  /// @brief Adds hits on surfaces and stores them
  /// @param hits vector of hits on surfaces
  /// @note The structure of @p hits is hits[Surfaces][Hits on a surface]
  void
  addHits(std::vector<std::vector<PlanarModuleCluster const*>>& hits) override;

  /// @brief Adds a hit structure to the list of hits
  /// @note This function does not test what is stored in the new element
  /// @param hit element added to the list
  void
  addSpacePoint(SpacePoint& sPoint) override;

  /// @brief Calculates the space points out of a given collection of hits and
  /// stores the data
  void
  calculateSpacePoints() override;

  /// @brief Returns the list of all stored space points
  /// @note This function is a pure getter of the current state
  /// @return full collection of all stored space points
  const std::vector<SpacePoint>&
  spacePoints() override;

protected:
  /// @brief Getter method for the local coordinates of a hit
  /// on its corresponding surface
  /// @param hit object related to the hit that holds the necessary information
  /// @return vector of the local coordinates of the hit on the surface
  Vector2D
  localCoords(const PlanarModuleCluster& hit) const override;

  /// @brief Getter method for the global coordinates of a hit
  /// @param hit object related to the hit that holds the necessary information
  /// @return vector of the global coordinates of the hit
  Vector3D
  globalCoords(const PlanarModuleCluster& hit) const override;

private:
  /// Storage of all stored data
  std::vector<SpacePoint> m_allSpacePoints;
};

}  // namespace Acts

#endif  // ACTS_TOOLS_ONEHITSPACEPOINTFINDER_H
