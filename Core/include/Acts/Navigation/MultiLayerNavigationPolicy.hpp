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

namespace Acts {

/// Concept for the indexed grid type
template <typename T>
concept IndexGridConcept = requires(T igrid) {
  typename T::grid_type;
  // grid must be of that type
  { igrid.grid } -> std::same_as<typename T::grid_type&>;
  // transform must be Transform3
  { igrid.transform } -> std::same_as<Transform3&>;
};

/// A navigation policy that uses grid based navigation for indexed surfaces
/// Navigate through a multilayer structure by creating an artificial path on
/// the grid.
template <IndexGridConcept indexed_grid>
class MultiLayerNavigationPolicy : public INavigationPolicy {
 public:
  struct Config {
    // The binning expansion for grid neighbor lookups
    std::vector<std::size_t> binExpansion = {0u, 0u};

    // The proto axis for the grid creation
    std::vector<DirectedProtoAxis> axis = {};
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
                                      const Logger& logger, Config config,
                                      indexed_grid grid)
      : m_volume(volume), m_indexedGrid(grid) {
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

    Acts::Experimental::detail::CenterReferenceGenerator rGenerator;
    Acts::Experimental::detail::IndexedGridFiller filler{config.binExpansion};
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
    const Transform3& transform = m_volume.transform();
    const Acts::Vector3 locPosition = transform.inverse() * args.position;
    const Acts::Vector3 locDirection = transform.linear() * args.direction;

    std::vector<Vector2> path = generatePath(locPosition, locDirection);

    const auto& surfaces = m_volume.surfaces();
    std::vector<const Surface*> surfCandidates = {};
    surfCandidates.reserve(surfaces.size());

    for (const auto& pos : path) {
      std::vector<std::size_t> indices = m_indexedGrid.grid.atPosition(pos);

      std::ranges::transform(indices, std::back_inserter(surfCandidates),
                             [&](const auto& i) { return &surfaces[i]; });
    }

    // remove duplicated candidates
    resolveDuplicates(surfCandidates);
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

    for (std::size_t i = 0; i < maxYIndex; i++) {
      auto v1 = m_indexedGrid.grid.lowerLeftBinEdge({1, i + 1});
      auto v2 = m_indexedGrid.grid.upperRightBinEdge({maxXIndex, i + 1});

      auto intersection = detail::IntersectionHelper2D::intersectSegment(
          Vector2(v1[0], v1[1]), Vector2(v2[0], v2[1]),
          Vector2(startPosition.x(), startPosition.y()),
          Vector2(direction.x(), direction.y()));
      if (!intersection.isValid()) {
        continue;
      }

      path.push_back(intersection.position());
    }
    return path;
  }
  /// Resolve duplicate on surface candidates
  /// @param surfaces is the surface candidates to check and resolve for duplicates
  void resolveDuplicates(std::vector<const Acts::Surface*>& surfaces) const {
    // sorting the surfaces according to their memory address
    std::ranges::sort(surfaces,
                      [](const Surface* a, const Surface* b) { return a < b; });

    surfaces.erase(std::unique(surfaces.begin(), surfaces.end()),
                   surfaces.end());
  }

 private:
  // The tracking volume
  const TrackingVolume& m_volume;

  // The grid that holds the indexed surfaces
  indexed_grid m_indexedGrid;
};

// Make alias for the static assert check for the Navigation Policy
using gridType =
    Grid<std::vector<std::size_t>,
         Axis<AxisType::Equidistant, Acts::AxisBoundaryType::Bound>,
         Axis<AxisType::Equidistant, Acts::AxisBoundaryType::Bound>>;
using indexedUpdatorType = Experimental::IndexedSurfacesNavigation<gridType>;

static_assert(
    NavigationPolicyConcept<MultiLayerNavigationPolicy<indexedUpdatorType>>);

}  // namespace Acts
