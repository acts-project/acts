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
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <array>
#include <memory>

namespace Acts {
namespace Experimental {

/// Helper struct that allows to fill surfaces into the candidate vector it
/// allows to use common navigation structs for volume, portal, surfaces
///
/// @param nState the navigation state
/// @param surfaces the surfaces that are filled in
inline static void fillSurfaceCandidates(
    const NavigationState& nState, const std::vector<const Surface*>& surfaces,
    NavigationState::SurfaceCandidates& candidates) {
  for (const auto& surface : surfaces) {
    candidates.push_back(NavigationState::SurfaceCandidate{
        ObjectIntersection<Surface>{}, surface, nullptr,
        nState.surfaceBoundaryCheck});
  }
}

/// Helper struct that allows to fill surfaces into the candidate vector it
/// allows to use common navigation structs for volume, portal, surfaces
///
/// @param nState the navigation state
/// @param portals the portals that are filled in
inline static void fillSurfaceCandidates(
    const NavigationState& /*nState*/,
    const std::vector<const Portal*>& portals,
    NavigationState::SurfaceCandidates& candidates) {
  for (const auto& portal : portals) {
    candidates.push_back(NavigationState::SurfaceCandidate{
        ObjectIntersection<Surface>{}, nullptr, portal, true});
  }
}

template <typename... updators_t>
class ChainedSurfaceCandidatesDelegate final
    : public ISurfaceCandidatesDelegate {
 public:
  /// The stored updators
  std::tuple<updators_t...> updators;

  /// Constructor for chained updators in a tuple, this will unroll
  /// the tuple and call them in sequence
  ///
  /// @param upts the updators to be called in chain
  ChainedSurfaceCandidatesDelegate(std::tuple<updators_t...> upts)
      : updators(std::move(upts)) {}

  /// A combined navigation state updator w/o intersection specifics
  ///
  /// @param gctx is the Geometry context of this call
  /// @param nState the navigation state to which the objects are attached
  void update(const GeometryContext& gctx, const NavigationState& nState,
              NavigationState::SurfaceCandidates& candidates) const final {
    // Unfold the tuple and add the attachers
    std::apply(
        [&](auto&&... updator) {
          (updator.update(gctx, nState, candidates), ...);
        },
        updators);
  }
};

}  // namespace Experimental
}  // namespace Acts
