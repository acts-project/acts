// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <unordered_map>
#include <vector>

namespace Acts {
namespace detail {

/// Clusterise tracks based on shared hits
///
/// @param trackMap : Multimap storing pair of track ID and vector of measurement ID. The keys are the number of measurement and are just there to focilitate the ordering.
/// @return an unordered map representing the clusters, the keys the ID of the primary track of each cluster and the store a vector of track IDs.
std::unordered_map<int, std::vector<int>> clusterDuplicateTracks(
    const std::multimap<int, std::pair<int, std::vector<int>>>& trackMap);

}  // namespace detail
}  // namespace Acts
