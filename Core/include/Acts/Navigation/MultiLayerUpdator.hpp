// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Navigation/NavigationStateUpdators.hpp"

#include <array>
#include <memory>

namespace Acts {
namespace Experimental {

template <typename grid_type, typename path_generator>
class MultiLayerSurfacesUpdatorImpl : public INavigationDelegate {
 public:
  /// The grid where the indices are stored
  grid_type grid;

  /// These are the cast parameters - copied from constructor
  std::array<BinningValue, grid_type::DIM> casts{};

  /// The updator, it can be indexed surfaces updator or indexed volumes updator
  updator_type updator;

  /// @brief  Constructor for a grid based surface attacher
  ///@param igrid the grid that is moved into this attacher
  /// @param icasts is the cast values array
  /// @param itr a transform applied to the global position
  MultiLayerSurfacesUpdatorImpl(
      grid_type&& igrid, const std::array<BinningValue, grid_type::DIM>& icasts,
      const Transform3& itr = Transform3::Identity())
      : grid(std::move(igrid)), casts(icasts), transform(itr) {}

  MultiLayerSurfacesUpdatorImpl() = delete;

  void update(const GeometryContext& gctx, NavigationState& nState) {
    auto path = path_generator(grid, nState);

    for (const auto& p : path) {
      const auto& entry = grid.atPosition(castPosition(p));
      auto extracted = IndexedSurfacesExtractor::extract(gctx, nState, entry);
    SurfacesFiller:
      fill(nState, extracted);
    }
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

}  // namespace Experimental

}  // namespace Acts
