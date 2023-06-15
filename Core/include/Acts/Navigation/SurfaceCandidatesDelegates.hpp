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
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/detail/DetectorVolumeExtractors.hpp"
#include "Acts/Detector/detail/IndexedLookupHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationDelegateHelpers.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>
#include <tuple>

namespace Acts {
namespace Experimental {

struct AllPortals final : public ISurfaceCandidatesDelegate {
  /// A ordered portal provider
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState is the navigation state to be updated
  ///
  /// @note that the intersections are ordered, such that the
  /// smallest intersection pathlength >= overstep tolerance is the lowest
  ///
  /// @return an ordered list of portal candidates
  inline void update(
      const GeometryContext& /*gctx*/, const NavigationState& nState,
      NavigationState::SurfaceCandidates& candidates) const final {
    auto currentVolume = nState.currentVolume;

    if (currentVolume == nullptr) {
      throw std::runtime_error("no detector volume set");
    }

    // Fill child volume portals
    for (const auto v : currentVolume->volumes()) {
      fillSurfaceCandidates(nState, v->portals(), candidates);
    }

    // Fill portals from the current volume
    fillSurfaceCandidates(nState, currentVolume->portals(), candidates);
  }
};

struct AllPortalsAndSurfaces final : public ISurfaceCandidatesDelegate {
  /// An ordered list of portals and surfaces provider
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState is the navigation state to be updated
  ///
  /// @note that the intersections are ordered, such that the
  /// smallest intersection pathlength >= overstep tolerance is the lowest
  ///
  /// @return an ordered list of portal and surface candidates
  inline void update(
      const GeometryContext& /*gctx*/, const NavigationState& nState,
      NavigationState::SurfaceCandidates& candidates) const final {
    auto currentVolume = nState.currentVolume;

    if (currentVolume == nullptr) {
      throw std::runtime_error("no detector volume set");
    }

    // Fill child volume portals
    for (const auto v : currentVolume->volumes()) {
      fillSurfaceCandidates(nState, v->portals(), candidates);
    }

    // Fill surfaces and portals from the current volume
    fillSurfaceCandidates(nState, currentVolume->portals(), candidates);
    fillSurfaceCandidates(nState, currentVolume->surfaces(), candidates);
  }
};

/// @brief This holds and extracts a collection of surfaces without much
/// checking, this could be e.g. support surfaces for layer structures,
/// e.g.
///
struct AdditionalSurfaces final : public ISurfaceCandidatesDelegate {
  /// The volumes held by this collection
  std::vector<const Surface*> surfaces = {};

  /// Extract the sub volumes from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call (ignored)
  /// @param nState is the navigation state
  ///
  inline void update(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState,
      NavigationState::SurfaceCandidates& candidates) const final {
    fillSurfaceCandidates(nState, surfaces, candidates);
  }
};

template <typename indexed_lookup_t>
struct GenericIndexedSurfaces : public ISurfaceCandidatesDelegate {
  indexed_lookup_t indexedLookup;

  GenericIndexedSurfaces(indexed_lookup_t _indexedLookup)
      : indexedLookup(std::move(_indexedLookup)) {}

  void update(const GeometryContext& gctx, const NavigationState& nState,
              NavigationState::SurfaceCandidates& candidates) const final {
    // Extract the index grid entry
    auto&& lookup = indexedLookup.lookup(gctx, nState);
    fillSurfaceCandidates(nState, lookup, candidates);
  }
};

/// @brief An indexed surface implementation access
///
/// @tparam grid_type is the grid type used for this indexed lookup
template <typename grid_type>
using IndexedSurfacesLookupHelper =
    IndexedLookupHelper<grid_type, IndexedSurfacesExtractor, SurfacesFiller>;

template <typename grid_type>
using IndexedSurfaces =
    GenericIndexedSurfaces<IndexedSurfacesLookupHelper<grid_type>>;

/// @brief An indexed surface implementaion with portal access
template <typename grid_type>
using IndexedSurfacesAllPortals =
    ChainedSurfaceCandidatesDelegate<IndexedSurfaces<grid_type>, AllPortals>;

template <typename grid_type>
inline static IndexedSurfaces<grid_type> makeIndexedSurfaces(
    grid_type igrid, const std::array<BinningValue, grid_type::DIM>& icasts,
    const Transform3& itr = Transform3::Identity()) {
  auto lookupHelper = IndexedSurfacesLookupHelper<grid_type>(
      std::forward<grid_type>(igrid), icasts, itr);
  return IndexedSurfaces<grid_type>(std::move(lookupHelper));
}

}  // namespace Experimental
}  // namespace Acts
