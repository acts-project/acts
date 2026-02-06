// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"

#include "Acts/Geometry/ReferenceGenerators.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"

namespace Acts::Experimental {

MultiLayerNavigationPolicy::MultiLayerNavigationPolicy(
    const GeometryContext& gctx, const TrackingVolume& volume,
    const Logger& logger, const Config& config, IndexedUpdatorType grid)
    : m_volume(volume), m_indexedGrid(std::move(grid)) {
  ACTS_VERBOSE("Constructing MultiLayerNavigationPolicy for volume "
               << m_volume.volumeName());

  // Fill the grid with surfaces
  std::vector<std::shared_ptr<const Surface>> surfaces = {};
  for (const auto& surface : m_volume.surfaces()) {
    if (!surface.isSensitive()) {
      continue;
    }
    surfaces.push_back(surface.getSharedPtr());
  }

  CenterReferenceGenerator rGenerator;
  IndexGridFiller filler{config.binExpansion};
  filler.fill(gctx, m_indexedGrid, surfaces, rGenerator, {});
}

void MultiLayerNavigationPolicy::initializeCandidates(
    const GeometryContext& gctx, const NavigationArguments& args,
    AppendOnlyNavigationStream& stream, const Logger& logger) const {
  ACTS_VERBOSE("MultiLayerNavigationPolicy Candidates initialization for volume"
               << m_volume.volumeName());
  const Transform3& itransform = m_volume.globalToLocalTransform(gctx);
  const Vector3 locPosition = itransform * args.position;
  const Vector3 locDirection = itransform.linear() * args.direction;

  std::vector<Vector2> path = generatePath(locPosition, locDirection);

  const auto& surfaces = m_volume.surfaces();
  std::vector<const Surface*> surfCandidates = {};
  surfCandidates.reserve(surfaces.size());

  for (const auto& pos : path) {
    // Local access
    std::vector<std::size_t> fAccessor = {0u, 1u};
    const auto& indices = m_indexedGrid.grid.atPosition(
        GridAccessHelpers::accessLocal<GridType>(pos, fAccessor));
    std::ranges::transform(indices, std::back_inserter(surfCandidates),
                           [&](const auto& i) { return &surfaces[i]; });
  }

  ACTS_VERBOSE("MultiLayerNavigationPolicy Candidates reported "
               << surfCandidates.size() << " candidates");

  // fill the navigation stream with the container
  for (const auto* surf : surfCandidates) {
    stream.addSurfaceCandidate(*surf, args.tolerance);
  }
}

std::vector<Vector2> MultiLayerNavigationPolicy::generatePath(
    const Vector3& startPosition, const Vector3& direction) const {
  std::vector<Vector2> path;

  auto maxXIndex = m_indexedGrid.grid.numLocalBins()[0];
  auto maxYIndex = m_indexedGrid.grid.numLocalBins()[1];
  Vector3 unitDir = direction.normalized();

  // cast the starting position and direction to the correct axis
  Vector2 startPoint{
      VectorHelpers::cast(startPosition, m_indexedGrid.casts[0]),
      VectorHelpers::cast(startPosition, m_indexedGrid.casts[1])};
  Vector2 startDir{VectorHelpers::cast(unitDir, m_indexedGrid.casts[0]),
                   VectorHelpers::cast(unitDir, m_indexedGrid.casts[1])};

  for (std::size_t i = 0; i < maxYIndex; i++) {
    auto v1 = m_indexedGrid.grid.lowerLeftBinEdge({1, i + 1});
    auto v2 = m_indexedGrid.grid.upperRightBinEdge({maxXIndex, i + 1});

    auto intersection = Acts::detail::IntersectionHelper2D::intersectSegment(
        Vector2(v1[0], v1[1]), Vector2(v2[0], v2[1]), startPoint, startDir);
    if (!intersection.isValid()) {
      continue;
    }

    path.push_back(intersection.position());
  }
  return path;
}

}  // namespace Acts::Experimental
