// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Detector/detail/IndexedGridFiller.hpp"
#include "Acts/Detector/detail/ReferenceGenerators.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Utilities/Grid.hpp"

#pragma once

namespace Acts {

/// Concept for the indexed grid type
template <typename T>
concept IndexGridConcept = requires(T igrid) {
  igrid.grid;
  igrid.transform;
  igrid.extractor;
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

  /// Main constructor, which creates the indexed surfaces assigned to a grid
  /// that is also constructed here
  /// @note Expects that the grid is defined but not filled - it will be filled here with the surfaces assigned to the @p volume
  /// @param gctx The geometrycontext object
  /// @param volume The tracking volume holding the surfaces that will be the indexed objects
  /// @param config The configuration of the Navigation Policy
  /// @param logger A logging instance
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

    // Define step size and number of steps
    m_step = std::sqrt(std::pow(m_indexedGrid.grid.binWidth()[0], 2) +
                       std::pow(m_indexedGrid.grid.binWidth()[1], 2));
    m_nSteps = m_indexedGrid.grid.numLocalBins()[1];
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
    const Acts::Vector3 locPosition = transform * args.position;
    const Acts::Vector3 locDirection = transform.linear() * args.direction;

    std::vector<Vector3> path =
        generatePath(locPosition, locDirection, m_step, m_nSteps);

    const auto& surfaces = m_volume.surfaces();
    std::vector<const Surface*> surfCandidates = {};

    for (const auto& pos : path) {
      std::vector<const Surface*> eSurfaces;
      std::vector<std::size_t> indices = m_indexedGrid.grid.atPosition(pos);
      eSurfaces.reserve(indices.size());

      std::for_each(indices.begin(), indices.end(),
                    [&](const auto& i) { eSurfaces.push_back(&surfaces[i]); });

      surfCandidates.insert(surfCandidates.end(), eSurfaces.begin(),
                            eSurfaces.end());
    }

    // remove duplicated candidates using geometryId
    resolveDuplicates(surfCandidates);
    ACTS_VERBOSE("MultiLayerNavigationPolicy Candidates reported"
                 << surfCandidates.size() << " candidates");

    // fill tne avigation stream with the container
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
  /// @param stepSize The step size for the path
  /// @param numberOfSteps The number of steps to take
  /// @return A vector of positions along the path
  std::vector<Vector3> generatePath(const Vector3& startPosition,
                                    const Vector3& direction, double stepSize,
                                    std::size_t numberOfSteps) const {
    std::vector<Vector3> path;
    path.reserve(numberOfSteps);
    for (std::size_t i = 0; i < numberOfSteps; ++i) {
      path.push_back(startPosition + i * stepSize * direction);
    }
    return path;
  }
  /// Resolve duplicate on surface candidates
  ///
  /// @param gctx The geometry context of the current geometry
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

  // The step size in the grid
  double m_step;

  // The number of steps in the grid
  std::size_t m_nSteps;
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
