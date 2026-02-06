// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/IndexGrid.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Surfaces/detail/IntersectionHelper2D.hpp"
#include "Acts/Utilities/Grid.hpp"

namespace Acts::Experimental {

/// A navigation policy that uses grid based navigation for indexed surfaces
/// Navigate through a multilayer structure by creating an artificial path on
/// the grid.
class MultiLayerNavigationPolicy : public INavigationPolicy {
 public:
  /// Type alias for 2D equidistant grid holding surface indices
  using GridType = Grid<std::vector<std::size_t>,
                        Axis<AxisType::Equidistant, AxisBoundaryType::Bound>,
                        Axis<AxisType::Equidistant, AxisBoundaryType::Bound>>;

  /// Type alias for indexed surfaces navigation updater
  using IndexedUpdatorType = IndexGrid<GridType>;

  /// Configuration for multilayer navigation behavior.
  struct Config {
    /// The binning expansion for grid neighbor lookups
    std::vector<std::size_t> binExpansion = {0u, 0u};
  };

  /// Main constructor, which expects the grid and will fill it with the
  /// surfaces from the volume passed
  /// @note Expects that the grid is defined but not filled - it will be filled here with the surfaces assigned to the @p volume
  /// @param gctx The geometrycontext object
  /// @param volume The tracking volume holding the surfaces that will be the indexed objects
  /// @param config The configuration of the Navigation Policy
  /// @param logger A logging instance
  /// @param grid The grid that will be filled with the surfaces
  explicit MultiLayerNavigationPolicy(const GeometryContext& gctx,
                                      const TrackingVolume& volume,
                                      const Logger& logger,
                                      const Config& config,
                                      IndexedUpdatorType grid);

  /// Update the navigation state from the surface array
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param stream The navigation stream to update
  /// @param logger The logger
  void initializeCandidates(const GeometryContext& gctx,
                            const NavigationArguments& args,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const;

  /// Connect this policy with a navigation delegate
  /// @param delegate The navigation delegate to connect to
  void connect(NavigationDelegate& delegate) const override {
    connectDefault<MultiLayerNavigationPolicy>(delegate);
  }

  /// Generate a path in the multilayer
  /// @param startPosition The starting position of the path (in local frame)
  /// @param direction The direction of the path (in local frame)
  /// @return A vector of positions along the path
  std::vector<Vector2> generatePath(const Vector3& startPosition,
                                    const Vector3& direction) const;

 private:
  // The tracking volume
  const TrackingVolume& m_volume;

  // The grid that holds the indexed surfaces
  IndexedUpdatorType m_indexedGrid;
};

static_assert(NavigationPolicyConcept<MultiLayerNavigationPolicy>);

}  // namespace Acts::Experimental
