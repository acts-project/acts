// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/Conversion/SeedConversion.hpp"

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

SimSeedContainer convertTracccToActsSeeds(
    std::vector<traccc::seed, std::pmr::polymorphic_allocator<traccc::seed>>&
        seeds,
    const std::map<std::size_t, std::size_t>& tracccToActsSpacepointIndexMap,
    const std::vector<SimSpacePoint>& actsSpacepoints) {
  auto fn = [&tracccToActsSpacepointIndexMap, actsSpacepoints](auto& seed) {
    return SimSeed(
        actsSpacepoints.at(tracccToActsSpacepointIndexMap.at(seed.spB_link)),
        actsSpacepoints.at(tracccToActsSpacepointIndexMap.at(seed.spM_link)),
        actsSpacepoints.at(tracccToActsSpacepointIndexMap.at(seed.spT_link)));
  };

  SimSeedContainer outputContainer;

  std::ranges::transform(seeds, std::back_inserter(outputContainer), fn);

  return outputContainer;
}

std::vector<traccc::seed, std::pmr::polymorphic_allocator<traccc::seed>>
convertActsToTracccSeeds(
    const SimSeedContainer& seeds,
    const SimSpacePointContainer& actsSpacepoints,
    const std::map<std::size_t, std::size_t>& actsToTracccSpacepointIndexMap) {
  auto fn = [&actsToTracccSpacepointIndexMap,
             &actsSpacepoints](const SimSeed& seed) {
    return traccc::seed{
        actsToTracccSpacepointIndexMap.at(
            std::distance(actsSpacepoints.data(), seed.sp()[0])),
        actsToTracccSpacepointIndexMap.at(
            std::distance(actsSpacepoints.data(), seed.sp()[0])),
        actsToTracccSpacepointIndexMap.at(
            std::distance(actsSpacepoints.data(), seed.sp()[0])),
        static_cast<typename traccc::point3::value_type>(seed.z()),
        seed.seedQuality()};
  };

  std::vector<traccc::seed, std::pmr::polymorphic_allocator<traccc::seed>>
      outputContainer;

  std::ranges::transform(seeds, std::back_inserter(outputContainer), fn);

  return outputContainer;
}

}  // namespace ActsExamples::Traccc::Common::Conversion
