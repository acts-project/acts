// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/IndexGrid.hpp"
#include "Acts/Geometry/ReferenceGenerators.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

namespace Acts {

class IndexGridNavigationConfig {
 public:
  /// The binning expansion for grid neighbor lookups
  std::vector<std::size_t> binExpansion = {};

  /// The reference expansion
  std::vector<double> referenceExpansion = {};

  /// A potential lookup surface - this would be intersected first,
  /// assumption is always closest forward
  std::shared_ptr<Surface> surface = nullptr;

  /// The reference generator
  std::shared_ptr<IReferenceGenerator> referenceGenerator =
      std::make_shared<PolyhedronReferenceGenerator>();
};

/// A navigation policy that uses grid based navigation for indexed surfaces
/// Navigate through a multilayer structure by creating an artificial path on
/// the grid.
template <typename GridType>
class IndexGridNavigationPolicy : public INavigationPolicy {
 public:
  using IndexGridType = IndexGrid<GridType>;

  /// Main constructor, which expects the grid and will fill it with the
  /// surfaces from the volume passed
  /// @note Expects that the grid is defined but not filled - it will be filled here with the surfaces assigned to the @p volume
  /// @param gctx The geometrycontext object
  /// @param volume The tracking volume holding the surfaces that will be the indexed objects
  /// @param logger A logging instance
  /// @param config The configuration of the Navigation Policy
  /// @param grid The index grid to use for navigation
  explicit IndexGridNavigationPolicy(const GeometryContext& gctx,
                                     const TrackingVolume& volume,
                                     const Logger& logger,
                                     const IndexGridNavigationConfig& config,
                                     const IndexGridType& grid)
      : m_cfg(config), m_volume(volume), m_indexGrid(grid) {
    ACTS_VERBOSE("Constructing IndexGridNavigationPolicy for volume '"
                 << m_volume.volumeName() << "'");

    // Fill the grid with the surfaces from the volume
    IndexGridFiller filler{m_cfg.binExpansion, m_cfg.referenceExpansion,
                           getDefaultLogger("IndexGridFiller", logger.level())};
    const auto& surfaces = m_volume.surfaces();
    // Fill the grid with surfaces
    std::vector<std::shared_ptr<const Surface>> surfacePtrs = {};
    surfacePtrs.reserve(surfaces.size());
    for (const auto& surface : surfaces) {
      surfacePtrs.push_back(surface.getSharedPtr());
    }
    filler.fill(gctx, m_indexGrid, surfacePtrs, *m_cfg.referenceGenerator);
  }

  /// Update the navigation state from the surface array
  /// @param gctx the geometry context
  /// @param args The navigation arguments
  /// @param stream The navigation stream to update
  /// @param logger The logger
  void initializeCandidates(const GeometryContext& gctx,
                            const NavigationArguments& args,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const {
    ACTS_VERBOSE(
        "IndexGridNavigationPolicy: candidates initialization for volume '"
        << m_volume.volumeName() << "'");

    // Nominal or intersected position
    Vector3 position = args.position;
    if (m_cfg.surface) {
      auto multiIntersection = m_cfg.surface->intersect(
          gctx, position, args.direction, BoundaryTolerance::Infinite());
      position = multiIntersection.closestForward().position();
    }

    // The indexing is guaranteed
    const auto& surfaces = m_volume.surfaces();
    const auto& indices =
        m_indexGrid.grid.atPosition(GridAccessHelpers::castPosition<GridType>(
            m_indexGrid.transform * position, m_indexGrid.casts));
    // Fill the navigation stream with the container
    for (const auto& idx : indices) {
      stream.addSurfaceCandidate(surfaces[idx], args.tolerance);
    }
  }

  /// Connect this policy with a navigation delegate
  /// @param delegate The navigation delegate to connect to
  void connect(NavigationDelegate& delegate) const override {
    connectDefault<IndexGridNavigationPolicy>(delegate);
  }

  /// @brief  Give const access to the index grid
  /// @return
  const IndexGridType& indexGrid() const { return m_indexGrid; }

 private:
  // Keep the config
  IndexGridNavigationConfig m_cfg;
  // The tracking volume
  const TrackingVolume& m_volume;
  // The index grid
  IndexGridType m_indexGrid;
};

// Regular cylinder is in phi and z
using RegularCylinderIndexGrid =
    Grid<std::vector<std::size_t>,
         Axis<AxisType::Equidistant, AxisBoundaryType::Closed>,
         Axis<AxisType::Equidistant, AxisBoundaryType::Bound>>;
using RegularCylinderIndexGridNavigationPolicy =
    IndexGridNavigationPolicy<RegularCylinderIndexGrid>;

static_assert(
    NavigationPolicyConcept<RegularCylinderIndexGridNavigationPolicy>);

// Regular ring in phi
using RegularRingIndexGrid =
    Grid<std::vector<std::size_t>,
         Axis<AxisType::Equidistant, AxisBoundaryType::Closed>>;
using RegularRingIndexGridNavigationPolicy =
    IndexGridNavigationPolicy<RegularRingIndexGrid>;

// Regular disc is in r and phi
using RegularDiscIndexGrid =
    Grid<std::vector<std::size_t>,
         Axis<AxisType::Equidistant, AxisBoundaryType::Bound>,
         Axis<AxisType::Equidistant, AxisBoundaryType::Closed>>;
using RegularDiscIndexGridNavigationPolicy =
    IndexGridNavigationPolicy<RegularDiscIndexGrid>;

static_assert(NavigationPolicyConcept<RegularDiscIndexGridNavigationPolicy>);

// Regular planar grid is in x and y
using RegularPlaneIndexGrid =
    Grid<std::vector<std::size_t>,
         Axis<AxisType::Equidistant, AxisBoundaryType::Bound>,
         Axis<AxisType::Equidistant, AxisBoundaryType::Bound>>;
using RegularPlaneIndexGridNavigationPolicy =
    IndexGridNavigationPolicy<RegularPlaneIndexGrid>;
static_assert(NavigationPolicyConcept<RegularPlaneIndexGridNavigationPolicy>);

}  // namespace Acts
