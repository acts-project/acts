// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Tools/ISpacePointBuilder.hpp"

namespace Acts {

/// @brief Structure for easier bookkeeping of hits.
struct SingleHitSpacePoint
{
  /// Storage of the hit(s) on a surface
  PlanarModuleCluster const* hitModule;
  /// Storage of a space point. Zero vector indicates unset point
  Vector3D spacePoint = {0., 0., 0.};
};

/// @class OneHitSpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the hits strip detectors need further treatment. This class takes
/// the digitized hits on a strip or pixel detector element and provides the
/// corresponding space point.
///
template <>
class SpacePointBuilder<SingleHitSpacePoint, void>
{
public:
  /// Default constructor
  SpacePointBuilder<SingleHitSpacePoint, void>() = delete;

  /// @brief Adds hits on surfaces and stores them
  /// @param spacePointStorage storage of hits and the therewith resulting space
  /// points
  /// @param hits vector of hits on the surface
  static void
  addHits(std::vector<SingleHitSpacePoint>& spacePointStorage,
          const std::vector<Acts::PlanarModuleCluster const*>& hits);

  /// @brief This function is intended to use hits on multiple surfaces that
  /// need to be combined in order to calculate a resulting space point. Since
  /// this is not needed for this class this function is deleted.
  static void
  addHits(std::vector<SingleHitSpacePoint>& spacePointStorage,
          const std::vector<Acts::PlanarModuleCluster const*>& hits1,
          const std::vector<Acts::PlanarModuleCluster const*>& hits2,
          const std::shared_ptr<void>                          cfg)
      = delete;

  /// @brief Calculates the space points out of a given collection of hits and
  /// stores the results
  /// @param spacePointStorage storage of the data
  static void
  calculateSpacePoints(std::vector<SingleHitSpacePoint>& spacePointStorage);

  /// @brief This function is intended to calculate space points out of given
  /// collection of hits using a specific configuration. Since this is not
  /// needed for this class this function is deleted.
  static void
  calculateSpacePoints(std::vector<SingleHitSpacePoint>& spacePointStorage,
                       const std::shared_ptr<void>       cfg)
      = delete;

protected:
  /// @brief Getter method for the local coordinates of a hit
  /// on its corresponding surface
  /// @param hit object related to the hit that holds the necessary information
  /// @return vector of the local coordinates of the hit on the surface
  static Vector2D
  localCoords(const PlanarModuleCluster& hit);

  /// @brief Getter method for the global coordinates of a hit
  /// @param hit object related to the hit that holds the necessary information
  /// @return vector of the global coordinates of the hit
  static Vector3D
  globalCoords(const PlanarModuleCluster& hit);
};

}  // namespace Acts
