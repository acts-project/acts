// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/detail/DetectorVolumeExtractors.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <array>

namespace Acts {
namespace Experimental {

/// @brief This is an index grid based navigation state updator, it uses
/// an extractor type and a filler type to handle the navigation state
///
/// @note a transform is applied `p3l = transform * p3` in order to allow
/// shifted, transformed grids
///
/// It can be used for volumes, surfaces at convenience
///
/// @tparam grid_t is the type of the grid
/// @tparam extractor_type is the helper to extract the object
/// @tparam filler_type is the helper to fill the object into the nState
template <typename grid_t, typename extractor_type, typename filler_type>
class IndexedLookupHelper {
 public:
  /// Broadcast the grid type
  using grid_type = grid_t;

  /// An extractor helper to get the object(s) from the volume
  extractor_type extractor;

  /// The grid where the indices are stored
  grid_type grid;

  /// These are the cast parameters - copied from constructor
  std::array<BinningValue, grid_type::DIM> casts{};

  /// @brief Constructor for a grid based surface attacher
  /// @param igrid the grid that is moved into this attacher
  /// @param icasts is the cast values array
  /// @param itr a transform applied to the global position
  IndexedLookupHelper(grid_type igrid,
                      const std::array<BinningValue, grid_type::DIM>& icasts,
                      const Transform3& itr = Transform3::Identity())
      : grid(std::move(igrid)), casts(icasts), transform(itr) {}

  /// @brief updates the navigation state with objects from the grid according
  /// to the filling type AFTER applying `p3loc = transform * p3`
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the surfaces are attached
  ///
  /// @note this is attaching objects without intersecting nor checking
  auto lookup(const GeometryContext& gctx,
              const NavigationState& nState) const {
    // Extract the index grid entry
    const auto& entry = grid.atPosition(castPosition(nState.position));
    return extractor.extract(gctx, nState, entry);
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
  /// A transform to be applied to the position
  Transform3 transform = Transform3::Identity();

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
