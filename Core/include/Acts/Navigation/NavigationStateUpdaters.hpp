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
#include "Acts/Detector/Portal.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <array>
#include <memory>

namespace Acts::Experimental {

/// Helper method to update the candidates (portals/surfaces),
/// this can be called for initial surface/portal estimation,
/// but also during the navigation to update the current list
/// of candidates.
///
/// @param gctx is the Geometry context of this call
/// @param nState [in,out] is the navigation state to be updated
///
/// @todo for surfaces skip the non-reached ones, while keep for portals
inline void updateCandidates(const GeometryContext& gctx,
                             NavigationState& nState) {
  const auto& position = nState.position;
  const auto& direction = nState.direction;

  NavigationState::SurfaceCandidates nextSurfaceCandidates;

  for (NavigationState::SurfaceCandidate c : nState.surfaceCandidates) {
    // Get the surface representation: either native surface of portal
    const Surface& sRep =
        c.surface != nullptr ? *c.surface : c.portal->surface();

    // Get the intersection @todo make a templated intersector
    // TODO surface tolerance
    auto sIntersection = sRep.intersect(gctx, position, direction,
                                        c.boundaryCheck, s_onSurfaceTolerance);
    for (auto& si : sIntersection.split()) {
      c.objectIntersection = si;
      if (c.objectIntersection &&
          c.objectIntersection.pathLength() > nState.overstepTolerance) {
        nextSurfaceCandidates.emplace_back(c);
      }
    }
  }

  nState.surfaceCandidates = std::move(nextSurfaceCandidates);
}

/// @brief This sets a single object, e.g. single surface or single volume
/// into the navigation state
///
/// @tparam navigation_type distinguishes between internal and external navigation
/// @tparam object_type the type of the object to be filled
/// @tparam filler_type is the helper to fill the object into nState
///
template <typename navigation_type, typename object_type, typename filler_type>
class SingleObjectNavigation : public navigation_type {
 public:
  /// Convenience constructor
  /// @param so the single object
  SingleObjectNavigation(const object_type* so) : m_object(so) {
    if (so == nullptr) {
      throw std::invalid_argument("SingleObjectNavigation: object is nullptr");
    }
  }

  /// @brief updates the navigation state with a single object that is filled in
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the surfaces are attached
  ///
  /// @note this is attaching objects without intersecting nor checking
  void update([[maybe_unused]] const GeometryContext& gctx,
              NavigationState& nState) const {
    filler_type::fill(nState, m_object);
  }

  /// Const Access to the object
  const object_type* object() const { return m_object; }

 private:
  // The object to be filled in
  const object_type* m_object = nullptr;
};

/// @brief This uses state less extractor and fillers to manipulate
/// the navigation state
///
/// @tparam navigation_type distinguishes between internal and external navigation
/// @tparam extractor_type the helper to extract the objects from
/// @tparam filler_type is the helper to fill the object into nState
///
template <typename navigation_type, typename extractor_type,
          typename filler_type>
class StaticAccessNavigation : public navigation_type {
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
/// @tparam navigation_type distinguishes between internal and external navigation
/// @tparam grid_t is the type of the grid
/// @tparam extractor_type is the helper to extract the object
/// @tparam filler_type is the helper to fill the object into the nState
///
template <typename navigation_type, typename grid_t, typename extractor_type,
          typename filler_type>
class IndexedGridNavigation : public navigation_type {
 public:
  /// Broadcast the grid type
  using grid_type = grid_t;

  /// An extractor helper to get the object(s) from the volume
  extractor_type extractor;

  /// The grid where the indices are stored
  grid_type grid;

  /// These are the cast parameters - copied from constructor
  std::array<BinningValue, grid_type::DIM> casts{};

  /// A transform to be applied to the position
  Transform3 transform = Transform3::Identity();

  /// @brief  Constructor for a grid based surface attacher
  /// @param igrid the grid that is moved into this attacher
  /// @param icasts is the cast values array
  /// @param itr a transform applied to the global position
  IndexedGridNavigation(grid_type&& igrid,
                        const std::array<BinningValue, grid_type::DIM>& icasts,
                        const Transform3& itr = Transform3::Identity())
      : grid(std::move(igrid)), casts(icasts), transform(itr) {}

  IndexedGridNavigation() = delete;

  /// @brief updates the navigation state with objects from the grid according
  /// to the filling type AFTER applying `p3loc = transform * p3`
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the surfaces are attached
  ///
  /// @note this is attaching objects without intersecting nor checking
  void update(const GeometryContext& gctx, NavigationState& nState) const {
    // Extract the index grid entry
    const auto& entry =
        grid.atPosition(GridAccessHelpers::castPosition<grid_type>(
            transform * nState.position, casts));
    auto extracted = extractor.extract(gctx, nState, entry);
    filler_type::fill(nState, extracted);

    // If the delegate type is of type IInternalNavigation
    if constexpr (std::is_base_of_v<IInternalNavigation, navigation_type>) {
      // Update the candidates
      updateCandidates(gctx, nState);
    }
  }
};

/// This is a chained extractor/filler implementation
/// Since there is no control whether it is a static or
/// payload extractor, these have to be provided by a tuple
///
/// @tparam navigation_type distinguishes between internal and external navigation
/// @tparam updators_t the updators that will be called in sequence
///
template <typename navigation_type, typename... updators_t>
class ChainedNavigation : public navigation_type {
 public:
  /// The stored updators
  std::tuple<updators_t...> updators;

  /// Constructor for chained updators in a tuple, this will unroll
  /// the tuple and call them in sequence
  ///
  /// @param upts the updators to be called in chain
  ChainedNavigation(const std::tuple<updators_t...>&& upts)
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

}  // namespace Acts::Experimental
