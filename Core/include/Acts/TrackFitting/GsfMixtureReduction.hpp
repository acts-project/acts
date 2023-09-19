// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"

namespace Acts {

/// Very simple mixture reduction method: Just removes the components with the
/// smallest weight until the required number of components is reached
///
/// @param cmpCache the component collection
/// @param maxCmpsAfterMerge the number of components we want to reach
/// @param surface the surface type on which the components are (unused in this method)
void reduceMixtureLargestWeights(
    std::vector<Acts::Experimental::GsfComponent> &cmpCache,
    std::size_t maxCmpsAfterMerge, const Surface &surface);

/// Greedy component reduction algorithm. Reduces the components with the
/// minimal symmetric KL-distance until the required number of components is
/// reached. Should yield an optimal reduction result at the cost of relatively
/// large computing cost.
///
/// @param cmpCache the component collection
/// @param maxCmpsAfterMerge the number of components we want to reach
/// @param surface the surface type on which the components are (unused in this method)
void reduceMixtureWithKLDistance(
    std::vector<Acts::Experimental::GsfComponent> &cmpCache,
    std::size_t maxCmpsAfterMerge, const Surface &surface);

/// More aggressive version of the KL-distance based mixture reducer. Tries to
/// merge more components in one pass without updating all distances. Should be
/// faster, but may merge components that have not the minimal distance.
///
/// @param cmpCache the component collection
/// @param maxCmpsAfterMerge the number of components we want to reach
/// @param surface the surface type on which the components are (unused in this method)
void reduceMixtureWithKLDistanceAggressive(
    std::vector<Acts::Experimental::GsfComponent> &cmpCache,
    std::size_t maxCmpsAfterMerge, const Surface &surface);

}  // namespace Acts
