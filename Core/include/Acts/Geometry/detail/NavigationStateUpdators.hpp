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
#include "Acts/Geometry/NavigationDelegates.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <array>
#include <memory>

namespace Acts {
namespace Experimental {
namespace detail {

/// @brief  This sets a single object, e.g. single surface or single volume
/// @tparam object_type the type of the object to be filled
/// @tparam filler_type is the helper to fill the object into nState
template <typename object_type, typename filler_type>
class SingleObjectImpl : public IDelegateImpl {
 public:
  // The object to be filled in
  const object_type* object = nullptr;

  /// Convenience constructor
  /// @param so the single object
  SingleObjectImpl(const object_type* so) : object(so) {}

  /// @brief updates the navigation state with a single object that is filled in
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the surfaces are attached
  ///
  /// @note this is attaching objects without intersecting nor checking
  void operator()([[maybe_unused]] const GeometryContext& gctx,
                  NavigationState& nState) const {
    filler_type::fill(nState, object);
  }
};

/// @brief This uses state less extractor and fillers to manipulate
/// the navigation state
///
/// @tparam extractor_type the helper to extract the objects from
/// @tparam filler_type is the helper to fill the object into nState
template <typename extractor_type, typename filler_type>
class StaticUpdatorImpl : public IDelegateImpl {
 public:
  /// @brief updates the navigation state with a single object that is filled in
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the surfaces are attached
  ///
  /// @note this is attaching objects without intersecting nor checking
  void operator()([[maybe_unused]] const GeometryContext& gctx,
                  NavigationState& nState) const {
    auto extracted = extractor_type::extract(gctx, *nState.currentVolume);
    filler_type::fill(nState, extracted);
  }
};

/// @brief  This is an index grid based navigation state updator, it uses
/// an extractor type and a filler type to handle the navigation state
///
/// It can be used for volumes, surfaces at convenience
///
/// @tparam grid_type is the type of the grid
/// @tparam extractor_type is the helper to extract the object
/// @tparam filler_type is the helper to fill the object into the nState
template <typename grid_type, typename extractor_type, typename filler_type>
class IndexedUpdatorImpl : public IDelegateImpl {
 public:
  /// The grid where the indices are stored
  grid_type grid;
  /// An extractor helper to get the object(s) from the volume
  extractor_type extractor;
  /// A filler helper to fill the selected objects into the naavigation state
  filler_type filler;
  /// These are the cast parameters
  std::array<BinningValue, grid_type::DIM> casts;

  /// @brief  Constructor for a grid based surface attacher
  /// @param igrid the grid that is moved into this attacher
  /// @param icasts is the cast values array
  IndexedUpdatorImpl(grid_type&& igrid,
                     const std::array<BinningValue, grid_type::DIM>& icasts)
      : grid(std::move(igrid)), casts(icasts) {}

  /// @brief updates the navigation state with objects from the grid according
  /// to the filling type
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the surfaces are attached
  ///
  /// @note this is attaching objects without intersecting nor checking
  void operator()([[maybe_unused]] const GeometryContext& gctx,
                  NavigationState& nState) const {
    // Extract the index grid entry
    const auto& entry = grid.atPosition(castPosition(nState.position));
    auto extracted =
        extractor_type::extract(gctx, *nState.currentVolume, entry);
    filler_type::fill(nState, extracted);
  }

  /// Cast into a lookup position
  ///
  /// @param position is the position of the update call
  std::array<ActsScalar, grid_type::DIM> castPosition(
      const Vector3& position) const {
    std::array<ActsScalar, grid_type::DIM> casted;
    fillCasts(position, casted,
              std::make_integer_sequence<std::size_t, grid_type::DIM>{});
    return casted;
  }

 private:
  /// Unroll cast loop
  /// @param position is the position of the update call
  /// @param a is the array to be filled
  template <typename Array, std::size_t... idx>
  void fillCasts(const Vector3& position, Array& a,
                 std::index_sequence<idx...>) const {
    ((a[idx] = VectorHelpers::cast(position, casts[idx])), ...);
  }
};

/// This is a chained extractor/filler implementation
/// Since there is no control whether it is a static or
/// payload extractor, these have to be provided by a tuple
///
/// @tparam updators_t the updators that will be called in sequence
template <typename... updators_t>
class ChainedUpdatorImpl : public IDelegateImpl {
 public:
  // The stored updators
  std::tuple<updators_t...> updators;

  /// Constructor for chained updators in a tuple, this will unroll
  /// the tuple and call them in sequence
  ///
  /// @param upts the updators to be called in chain
  ChainedUpdatorImpl(const std::tuple<updators_t...>& upts) : updators(upts) {}

  /// A combined navigation state updator w/o intersection specifics
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the objects are attached
  ///
  void operator()(const GeometryContext& gctx, NavigationState& nState) const {
    // Unfold the tuple and add the attachers
    std::apply([&](auto&&... updator) { ((updator(gctx, nState)), ...); },
               updators);
  }
};

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
