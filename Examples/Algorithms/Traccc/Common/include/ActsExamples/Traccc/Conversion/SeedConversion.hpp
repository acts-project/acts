// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Seed.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <cstdint>
#include <cstdlib>
#include <memory_resource>

#include "traccc/edm/seed.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "vecmem/containers/vector.hpp"

namespace ActsExamples::Traccc::Common::Conversion {

/// @brief Converts a collection of traccc seeds and appends the result to the given outputContainer.
/// @param seeds the traccc seeds
/// @param SpacePointConv the spacepoint ConversionData (traccc space point -> acts space point).
/// @param outputContainer the container to put the converted space points into (an empty container is expected).
/// @returns The seed ConversionData (traccc seed -> acts seed).
SimSeedContainer convertTracccToActsSeeds(
    std::vector<traccc::seed, std::pmr::polymorphic_allocator<traccc::seed>>&
        seeds,
    const std::map<std::size_t, std::size_t>& tracccToActsSpacepointIndexMap,
    const std::vector<SimSpacePoint>& actsSpacepoints);

/// @brief Converts a collection of seeds to traccc seeds and appends the result to the given outputContainer.
/// @param seeds theseeds
/// @param SpacePointConv the spacepoint ConversionData (acts space point -> traccc space point).
/// @param outputContainer the container to put the converted space points into (an empty container is expected).
/// @returns The seed ConversionData (acts seed -> traccc seed).
std::vector<traccc::seed, std::pmr::polymorphic_allocator<traccc::seed>>
convertActsToTracccSeeds(
    const SimSeedContainer& seeds,
    const SimSpacePointContainer& actsSpacepoints,
    const std::map<std::size_t, std::size_t>& actsToTracccSpacepointIndexMap);

}  // namespace ActsExamples::Traccc::Common::Conversion
