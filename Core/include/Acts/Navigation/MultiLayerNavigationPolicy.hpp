// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/detail/IndexedGridFiller.hpp"
#include "Acts/Detector/detail/ReferenceGenerators.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Surfaces/detail/IntersectionHelper2D.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

namespace Acts::Experimental {

/// A navigation policy that uses grid based navigation for indexed surfaces
/// Navigate through a multilayer structure by creating an artificial path on
/// the grid.
class MultiLayerNavigationPolicy : public INavigationPolicy {
 public:
  using gridType =
      Grid<std::vector<std::size_t>,
           Axis<AxisType::Equidistant, Acts::AxisBoundaryType::Bound>,
           Axis<AxisType::Equidistant, Acts::AxisBoundaryType::Bound>>;
  using indexedUpdatorType = IndexedSurfacesNavigation<gridType>;

  struct Config {
    // The binning expansion for grid neighbor lookups
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
                                      indexedUpdatorType grid)
      : m_volume(volume), m_indexedGrid(std::move(grid)) {
    ACTS_VERBOSE("Constructing MultiLayerNavigationPolicy for volume "
                 << m_volume.volumeName());

    // Fill the grid with surfaces
    std::vector<std::shared_ptr<const Surface>> surfaces = {};
    for (const auto& surface : m_volume.surfaces()) {
      if (surface.associatedDetectorElement() == nullptr) {
        continue;
      }
      surfaces.push_back(surface.getSharedPtr());
    }

    Experimental::detail::CenterReferenceGenerator rGenerator;
    Experimental::detail::IndexedGridFiller filler{config.binExpansion};
    filler.fill(gctx, m_indexedGrid, surfaces, rGenerator, {});
  }

  /// Update the navigation state from the surface array
  /// @param args The navigation arguments
  /// @param stream The navigation stream to update
  /// @param logger The logger
  void initializeCandidates(const NavigationArguments& args,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const {
    ACTS_VERBOSE(
        "MultiLayerNavigationPolicy Candidates initialization for volume"
        << m_volume.volumeName());
    const Transform3& itransform = m_volume.itransform();
    const Vector3 locPosition = itransform * args.position;
    const Vector3 locDirection = itransform * args.direction;

    std::vector<Vector2> path = generatePath(locPosition, locDirection);

    const auto& surfaces = m_volume.surfaces();
    std::vector<const Surface*> surfCandidates = {};
    surfCandidates.reserve(surfaces.size());

    for (const auto& pos : path) {
      std::vector<std::size_t> indices = m_indexedGrid.grid.atPosition(pos);

      std::ranges::transform(indices, std::back_inserter(surfCandidates),
                             [&](const auto& i) { return &surfaces[i]; });
    }

    ACTS_VERBOSE("MultiLayerNavigationPolicy Candidates reported"
                 << surfCandidates.size() << " candidates");

    // fill the navigation stream with the container
    for (const auto* surf : surfCandidates) {
      stream.addSurfaceCandidate(*surf, args.tolerance);
    }
  }

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
                                    const Vector3& direction) const {
    std::vector<Vector2> path;

    auto maxXIndex = m_indexedGrid.grid.numLocalBins()[0];
    auto maxYIndex = m_indexedGrid.grid.numLocalBins()[1];
    Vector3 unitDir = direction.normalized();

    for (std::size_t i = 0; i < maxYIndex; i++) {
      auto v1 = m_indexedGrid.grid.lowerLeftBinEdge({1, i + 1});
      auto v2 = m_indexedGrid.grid.upperRightBinEdge({maxXIndex, i + 1});

      auto intersection = Acts::detail::IntersectionHelper2D::intersectSegment(
          Vector2(v1[0], v1[1]), Vector2(v2[0], v2[1]),
          startPosition.template block<2, 1>(0, 0),
          unitDir.template block<2, 1>(0, 0));
      if (!intersection.isValid()) {
        continue;
      }

      path.push_back(intersection.position());
    }
    return path;
  }

 private:
  // The tracking volume
  const TrackingVolume& m_volume;

  // The grid that holds the indexed surfaces
  indexedUpdatorType m_indexedGrid;
};

static_assert(NavigationPolicyConcept<MultiLayerNavigationPolicy>);

}  // namespace Acts::Experimental
