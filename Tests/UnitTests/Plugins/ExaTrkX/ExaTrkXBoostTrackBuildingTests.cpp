// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/ExaTrkX/BoostTrackBuilding.hpp"

#include <algorithm>
#include <random>

BOOST_AUTO_TEST_CASE(test_track_building) {
    std::mt19937 gen{42};
    
    // Make some spacepoint IDs
    std::vector<int> spacepointIds(16);
    std::iota(spacepointIds.begin(), spacepointIds.end(), 0);
    
    // Build 4 tracks with 4 hits
    std::vector<std::vector<int>> refTracks;
    for(auto t=0ul; t < 4; ++t) {
        refTracks.emplace_back(spacepointIds.begin()+4*t, spacepointIds.begin()+4*(t+1));
    }
    
    // Make edges
    std::vector<int> edges;
    for(const auto &track : refTracks) {
        for(auto it=track.begin(); it != track.end()-1; ++it) {
            edges.push_back(*it);
            edges.push_back(*std::next(it));
        }
    }
    
    
}
