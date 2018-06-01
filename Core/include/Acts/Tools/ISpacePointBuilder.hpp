// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include "Acts/Digitization/PlanarModuleCluster.hpp"

namespace Acts {

/// @struct SpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the hits on some detector elements need further treatment. This struct
/// serves as default structure of the process to take the digitized clusters on
/// a
/// detector element and provide the corresponding space point. The empty
/// structe is used to forbid the usage of an arbitrary data type as template
/// parameter and enforces the implementation of explicit structures.
///
/// @note The choice of which kind of data should be treated in which way is
/// steered by the choice of the template parameter. This parameter represents a
/// structure that needs to store at least a cluster/multiple clusters and the
/// corresponding space point. The second template parameter allows the
/// application of different configurations for the application.
///
template <class S, class C = void>
struct SpacePointBuilder
{
};

namespace SPB {
  /// @brief Adds clusters to the list
  /// @param spacePointStorage storage of clusters and the therewith resulting
  /// space
  /// points
  /// @param clusters list of hits
  /// @note This function is intended to be used for the case that a single
  /// cluster
  /// (e.g. in a pixel detector module) results in a space point.
  template <class S>
  static void
  addClusters(std::vector<S>& spacePointStorage,
              const std::vector<Acts::PlanarModuleCluster const*>& clusters)
  {
    SpacePointBuilder<S, void>::addClusters(spacePointStorage, clusters);
  }

  /// @brief Adds clusters to the list
  /// @param spacePointStorage storage of clusters and the therewith resulting
  /// space
  /// points
  /// @param clusters1 list of clusters
  /// @param clusters2 list of clusters on another surface(s) than @p clusters1
  /// @param cfg optional configuration to steer the combinations of the
  /// elements of @p clusters1 and @p clusters2
  /// @note This function is intended to be used for the case that two clusters
  /// (e.g. in a double strip detector module) result in a space point.
  template <class S, class C>
  static void
  addClusters(std::vector<S>& spacePointStorage,
              const std::vector<Acts::PlanarModuleCluster const*>& clusters1,
              const std::vector<Acts::PlanarModuleCluster const*>& clusters2,
              const std::shared_ptr<C> cfg = nullptr)
  {
    SpacePointBuilder<S, C>::addClusters(
        spacePointStorage, clusters1, clusters2, cfg);
  }

  /// @brief Calculates the space points out of a given collection of clusters
  /// and
  /// stores the results
  /// @param spacePointStorage storage of the clusters and the corresponding
  /// space
  /// points
  /// @param cfg optional configuration to steer the calculation of space points
  template <class S>
  static void
  calculateSpacePoints(std::vector<S>& spacePointStorage)
  {
    SpacePointBuilder<S, void>::calculateSpacePoints(spacePointStorage);
  }

  /// @brief Calculates the space points out of a given collection of clusters
  /// and
  /// stores the results
  /// @param spacePointStorage storage of the clusters and the corresponding
  /// space
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
