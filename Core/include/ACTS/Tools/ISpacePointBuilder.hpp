// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include "ACTS/Digitization/PlanarModuleCluster.hpp"

namespace Acts {

/// @struct SpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the hits on some detector elements need further treatment. This struct
/// serves as default structure of the process to take the digitized hits on a
/// detector element and provide the corresponding space point. The empty
/// structe is used to forbid the usage of an arbitrary data type as template
/// parameter and enforces the implementation of explicit structures.
///
/// @note The choice of which kind of data should be treated in which way is
/// steered by the choice of the template parameter. This parameter represents a
/// structure that needs to store at least a hit/multiple hits and the
/// corresponding space point. The second template parameter allows the
/// application of different configurations for the application.
///
template <class S, class C = void>
struct SpacePointBuilder
{
};

namespace SPB {
  /// @brief Adds hits to the list
  /// @param spacePointStorage storage of hits and the therewith resulting space
  /// points
  /// @param hits list of hits
  /// @note This function is intended to be used for the case that a single hit
  /// (e.g. in a pixel detector module) results in a space point.
  template <class S>
  static void
  addHits(std::vector<S> spacePointStorage,
          const std::vector<Acts::PlanarModuleCluster const*>& hits)
  {
    SpacePointBuilder<S, void>::addHits(spacePointStorage, hits);
  }

  /// @brief Adds hits to the list
  /// @param spacePointStorage storage of hits and the therewith resulting space
  /// points
  /// @param hits1 list of hits
  /// @param hits2 list of hits on another surface(s) than @p hits1
  /// @param cfg optional configuration to steer the combinations of the
  /// elements of @p hits1 and @p hits2
  /// @note This function is intended to be used for the case that two hits
  /// (e.g. in a double strip detector module) result in a space point.
  template <class S, class C>
  static void
  addHits(std::vector<S> spacePointStorage,
          const std::vector<Acts::PlanarModuleCluster const*>& hits1,
          const std::vector<Acts::PlanarModuleCluster const*>& hits2,
          const std::shared_ptr<C>                             cfg = nullptr)
  {
    SpacePointBuilder<S, C>::addHits(spacePointStorage, hits1, hits2, cfg);
  }

  /// @brief Calculates the space points out of a given collection of hits and
  /// stores the results
  /// @param spacePointStorage storage of the hits and the corresponding space
  /// points
  /// @param cfg optional configuration to steer the calculation of space points
  template <class S>
  static void
  calculateSpacePoints(std::vector<S>& spacePointStorage)
  {
    SpacePointBuilder<S, void>::calculateSpacePoints(spacePointStorage);
  }

  /// @brief Calculates the space points out of a given collection of hits and
  /// stores the results
  /// @param spacePointStorage storage of the hits and the corresponding space
  /// points
  /// @param cfg optional configuration to steer the calculation of space points
  template <class S, class C>
  static void
  calculateSpacePoints(std::vector<S>&          spacePointStorage,
                       const std::shared_ptr<C> cfg = nullptr)
  {
    SpacePointBuilder<S, C>::calculateSpacePoints(spacePointStorage, cfg);
  }

}  // namespace SP

}  // namespace Acts
