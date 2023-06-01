// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <array>
#include <memory>

namespace Acts {
namespace Experimental {

/// @brief  This sets a single object, e.g. single surface or single volume
/// @tparam object_type the type of the object to be filled
/// @tparam filler_type is the helper to fill the object into nState
template <typename object_type, typename filler_type>
class SingleObjectImpl : public INavigationDelegate {
 public:
  /// Convenience constructor
  /// @param so the single object
  SingleObjectImpl(const object_type* so) : object(so) {}

  /// @brief updates the navigation state with a single object that is filled in
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the surfaces are attached
  ///
  /// @note this is attaching objects without intersecting nor checking
  void update([[maybe_unused]] const GeometryContext& gctx,
              NavigationState& nState) const {
    filler_type::fill(nState, object);
  }

 private:
  // The object to be filled in
  const object_type* object = nullptr;
};

/// @brief This uses state less extractor and fillers to manipulate
/// the navigation state
///
/// @tparam extractor_type the helper to extract the objects from
/// @tparam filler_type is the helper to fill the object into nState
template <typename extractor_type, typename filler_type>
class StaticUpdatorImpl : public INavigationDelegate {
 public:
  /// @brief updates the navigation state with a single object that is filled in
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the surfaces are attached
  ///
  /// @note this is attaching objects without intersecting nor checking
  void update([[maybe_unused]] const GeometryContext& gctx,
              NavigationState& nState) const {
    auto extracted = extractor_type::extract(gctx, nState);
    filler_type::fill(nState, extracted);
  }
};

/// @brief  This is an index grid based navigation state updator, it uses
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
class IndexedUpdatorImpl : public INavigationDelegate {
 public:
  /// Broadcast the grid type
  using grid_type = grid_t;

  /// An extractor helper to get the object(s) from the volume
  extractor_type extractor;

  /// The grid where the indices are stored
  grid_type grid;

  /// These are the cast parameters - copied from constructor
  std::array<BinningValue, grid_type::DIM> casts{};

  /// @brief  Constructor for a grid based surface attacher
  /// @param igrid the grid that is moved into this attacher
  /// @param icasts is the cast values array
  /// @param itr a transform applied to the global position
  IndexedUpdatorImpl(grid_type&& igrid,
                     const std::array<BinningValue, grid_type::DIM>& icasts,
                     const Transform3& itr = Transform3::Identity())
      : grid(std::move(igrid)), casts(icasts), transform(itr) {}

  IndexedUpdatorImpl() = delete;

  /// @brief updates the navigation state with objects from the grid according
  /// to the filling type AFTER applying `p3loc = transform * p3`
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the surfaces are attached
  ///
  /// @note this is attaching objects without intersecting nor checking
  void update(const GeometryContext& gctx, NavigationState& nState) const {
    // Extract the index grid entry
    const auto& entry = grid.atPosition(castPosition(nState.position));
    auto extracted = extractor.extract(gctx, nState, entry);
    filler_type::fill(nState, extracted);
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

/// This is a chained extractor/filler implementation
/// Since there is no control whether it is a static or
/// payload extractor, these have to be provided by a tuple
///
/// @tparam updators_t the updators that will be called in sequence
template <typename... updators_t>
class ChainedUpdatorImpl : public INavigationDelegate {
 public:
  /// The stored updators
  std::tuple<updators_t...> updators;

  /// Constructor for chained updators in a tuple, this will unroll
  /// the tuple and call them in sequence
  ///
  /// @param upts the updators to be called in chain
  ChainedUpdatorImpl(const std::tuple<updators_t...>&& upts)
      : updators(std::move(upts)) {}

  /// A combined navigation state updator w/o intersection specifics
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the objects are attached
  ///
  void update(const GeometryContext& gctx, NavigationState& nState) const {
    // Unfold the tuple and add the attachers
    std::apply(
        [&](auto&&... updator) { ((updator.update(gctx, nState)), ...); },
        updators);
  }
};

}  // namespace Experimental
}  // namespace Acts
