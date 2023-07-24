// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
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
/// @param surfaces the surfaces that are filled in
/// @param candidates the surface candidates to be updated
inline static void fillSurfaceCandidates(
    NavigationState::SurfaceCandidates& candidates,
    const std::vector<const Surface*>& surfaces, bool boundaryCheck) {
  for (const auto& surface : surfaces) {
    candidates.push_back(NavigationState::SurfaceCandidate{
        ObjectIntersection<Surface>{}, surface, nullptr, boundaryCheck});
  }
}

/// Helper struct that allows to fill surfaces into the candidate vector it
/// allows to use common navigation structs for volume, portal, surfaces
///
/// @param portals the portals that are filled in
/// @param candidates the surface candidates to be filled
inline static void fillSurfaceCandidates(
    NavigationState::SurfaceCandidates& candidates,
    const std::vector<const Portal*>& portals) {
  for (const auto& portal : portals) {
    candidates.push_back(NavigationState::SurfaceCandidate{
        ObjectIntersection<Surface>{}, nullptr, portal, true});
  }
}

}  // namespace Experimental
}  // namespace Acts
