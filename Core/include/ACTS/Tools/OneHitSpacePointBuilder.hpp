// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_TOOLS_ONEHITSPACEPOINTBUILDER_H
#define ACTS_TOOLS_ONEHITSPACEPOINTBUILDER_H

#include <ACTS/Digitization/PlanarModuleCluster.hpp>
#include <array>
#include "ACTS/Digitization/CartesianSegmentation.hpp"
#include "ACTS/Digitization/DigitizationModule.hpp"

namespace Acts {

/// @class OneHitSpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the hits strip detectors need further treatment. This class takes
/// the digitized hits and provides the corresponding space point.
///
class OneHitSpacePointBuilder
{
public:
  /// @brief Structure for easier bookkeeping of hits.
  struct Hit
  {
    /// Storage of the hit on a surface
    Acts::PlanarModuleCluster const* hitModule;
    /// Storage of a space point. Zero vector indicates unset point
    Vector3D spacePoint = {0., 0., 0.};
  };

  /// Constructor
  OneHitSpacePointBuilder();

  /// @brief Adds hits on surfaces and stores them
  /// @param vec vector of hits on a surface
  void
  addHits(const std::vector<Acts::PlanarModuleCluster>& vec);

  /// @brief Adds a hit structure to the list of hits
  /// @note This function does not test what is stored in the new element
  /// @param hit element added to the list
  void
  addHit(const Hit& hit);

  /// @brief Calculates the space points out of a given collection of hits and stores the data
  void
  calculateSpacePoints();

  /// @brief Returns the list of all stored space points
  /// @note This function is a pure getter of the current state
  /// @return full collection of all stored space points
  const std::vector<Hit>&
  hits();

private:
  /// Storage of all stored data
  std::vector<OneHitSpacePointBuilder::Hit> m_allHits;

  /// @brief Getter method for the local coordinates of a hit
  /// on its corresponding surface
  /// @param hit object related to the hit that holds the necessary information
  /// @return vector of the local coordinates of the hit on the surface
  Acts::Vector2D
  localCoords(const Acts::PlanarModuleCluster& hit) const;

  /// @brief Getter method for the global coordinates of a hit
  /// @param hit object related to the hit that holds the necessary information
  /// @return vector of the global coordinates of the hit
  Acts::Vector3D
  globalCoords(const Acts::PlanarModuleCluster& hit) const;
};

}  // namespace Acts

#endif  // ACTS_TOOLS_ONEHITSPACEPOINTFINDER_H
