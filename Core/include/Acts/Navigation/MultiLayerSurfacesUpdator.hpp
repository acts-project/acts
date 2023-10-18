// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Navigation/NavigationStateUpdators.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <array>
#include <memory>

namespace Acts {
namespace Experimental {

template <typename grid_t, typename path_generator>
class MultiLayerSurfacesUpdatorImpl : public INavigationDelegate {
 public:
  /// Broadcast the grid type
  using grid_type = grid_t;

  /// The grid where the indices are stored
  grid_type grid;

  /// These are the cast parameters - copied from constructor
  std::array<BinningValue, grid_type::DIM> casts{};

  /// A transform to be applied to the position
  Transform3 transform = Transform3::Identity();

  /// The path generator
  path_generator pgenerator;

  /// @brief  Constructor for a grid based surface attacher
  ///@param igrid the grid that is moved into this attacher
  /// @param icasts is the cast values array
  /// @param itr a transform applied to the global position
  MultiLayerSurfacesUpdatorImpl(
      grid_type&& igrid, const std::array<BinningValue, grid_type::DIM>& icasts,
      const Transform3& itr = Transform3::Identity())
      : grid(std::move(igrid)), casts(icasts), transform(itr) {}

  MultiLayerSurfacesUpdatorImpl() = delete;

  void update(const GeometryContext& gctx, NavigationState& nState) const {
    auto step = std::sqrt(std::pow(grid.binWidth()[0], 2) +
                          std::pow(grid.binWidth()[1], 2));
    auto path = pgenerator(nState.position, nState.direction, step,
                           grid.numLocalBins()[1]);

    std::vector<const Acts::Surface*> surfCandidates = {};

    for (const auto& p : path) {
      const auto& entry = grid.atPosition(castPosition(p));
      const auto extracted =
          IndexedSurfacesExtractor::extract(gctx, nState, entry);
      surfCandidates.insert(surfCandidates.end(), extracted.begin(),
                            extracted.end());
    }

    resolveDuplicates(gctx, surfCandidates);
    SurfacesFiller::fill(nState, surfCandidates);
  }

  /// Cast into a lookup position
  ///
  /// @param position is the position of the update call
  std::array<ActsScalar, grid_type::DIM> castPosition(
      const Vector3& position) const {
    // Transform into local 3D frame
    Vector3 tposition = transform * position;

    std::array<ActsScalar, grid_type::DIM> casted{};
    fillCasts(tposition, casted,
              std::make_integer_sequence<std::size_t, grid_type::DIM>{});
    return casted;
  }

  /// Resolve duplicate on surface candidates
  ///
  /// @param gctx The geometry context of the current geometry
  /// @param surfaces is the surface candidates to check and resolve for duplicates
  void resolveDuplicates(const GeometryContext& gctx,
                         std::vector<const Acts::Surface*>& surfaces) const {
    // sorting the surfaces according to their radial distance
    std::sort(surfaces.begin(), surfaces.end(),
              [&gctx](const auto& surf1, const auto& surf2) {
                if (surf1->center(gctx).x() != surf2->center(gctx).x()) {
                  return surf1->center(gctx).x() < surf2->center(gctx).x();
                }
                if (surf1->center(gctx).y() != surf2->center(gctx).y()) {
                  return surf1->center(gctx).y() < surf2->center(gctx).y();
                }
                return surf1->center(gctx).z() < surf2->center(gctx).z();
              });

    // Remove the duplicates
    surfaces.erase(std::unique(surfaces.begin(), surfaces.end()),
                   surfaces.end());
  }

 private:
  /// Unroll the cast loop
  /// @param position is the position of the update call
  /// @param a is the array to be filled
  template <typename Array, std::size_t... idx>
  void fillCasts(const Vector3& position, Array& a,
                 std::index_sequence<idx...> /*indices*/) const {
    ((a[idx] = VectorHelpers::cast(position, casts[idx])), ...);
  }
};

struct PathGridSurfacesGenerator {
  std::vector<Vector3> operator()(Vector3 startPosition,
                                  const Vector3& direction, ActsScalar stepSize,
                                  std::size_t numberOfSteps) const {
    std::vector<Vector3> pathCoordinates = {};
    pathCoordinates.reserve(numberOfSteps);

    auto tposition = std::move(startPosition);
    auto stepSizeY = stepSize * sin(Acts::VectorHelpers::phi(direction));
    auto stepSizeX = stepSize * cos(Acts::VectorHelpers::phi(direction));

    for (std::size_t i = 0; i < numberOfSteps; i++) {
      pathCoordinates.push_back(tposition);
      tposition.y() = tposition.y() + stepSizeY;
      tposition.x() = tposition.x() + stepSizeX;
    }

    return pathCoordinates;
  }
};

}  // namespace Experimental

}  // namespace Acts
