// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s)
#include "Acts/Seeding/Seed.hpp"

// Acts Examples include(s)
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Traccc/Common/Util/MapUtil.hpp"

// Traccc include(s)
#include "traccc/edm/seed.hpp"
#include "traccc/edm/spacepoint.hpp"

// VecMem include(s)
#include "vecmem/containers/vector.hpp"

// System include(s)
#include <cstdint>
#include <cstdlib>

namespace ActsExamples::Traccc::Common::Conversion {

template <typename map_t>
SimSeed convertSeed(const traccc::seed& seed, const traccc::spacepoint_collection_types::const_view& spacePointsView, const map_t& spacePointMap){
    const std::array<traccc::spacepoint, 3> spacePoints = seed.get_spacepoints(spacePointsView);
    return SimSeed(
        spacePointMap.at(spacePoints[0]),
        spacePointMap.at(spacePoints[1]),
        spacePointMap.at(spacePoints[2]),
        seed.z_vertex,
        seed.weight
    );
}

template <typename allocator_t, typename map_t>
auto convertSeeds(std::vector<traccc::seed, allocator_t>& seeds, const traccc::spacepoint_collection_types::const_view& spacePointsView, const map_t& spacePointMap){
    auto fn = [&spacePointsView, &spacePointMap](auto& seed){
        return convertSeed(seed, spacePointsView, spacePointMap);
    };
    return Util::convert<traccc::seed, SimSeed>(seeds, fn);
}
/*
template <typename map_t>
SimSeed convertSeed(SimSeed& seed, const map_t& spacePointMap){
    const std::array<traccc::spacepoint, 3> spacePoints = seed.get_spacepoints(spacePointsView);
    return traccc::seed{,,,seed.z(), seed.seedQuality()};
}

template <typename allocator_t, typename map_t>
auto convertSeeds(SimSeedContainer& seeds, const map_t& spacePointMap){
    auto fn = [&spacePointsView, &spacePointMap](auto& seed){
        return convertSeed(seed, spacePointsView, spacePointMap);
    };
    return Util::convert<traccc::seed, SimSeed>(seeds, fn);
}*/

}  // namespace Acts::TracccPlugin
