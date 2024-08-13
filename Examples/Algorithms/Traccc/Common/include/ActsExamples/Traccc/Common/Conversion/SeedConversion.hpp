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

// Custom hash and equals functions
// as some are not defined by std::hash and std::equal_to

struct TracccSeedHash{
    std::size_t operator()(const traccc::seed& s) const noexcept {
        return s.spB_link ^ s.spM_link ^ s.spT_link;
    }
};

struct TracccSeedEquals {
    bool operator()(const traccc::seed& lhs, const traccc::seed& rhs) const {
        return lhs.spB_link == rhs.spB_link &&
        lhs.spM_link == rhs.spM_link &&
        lhs.spT_link == rhs.spT_link &&
        lhs.weight == rhs.weight &&
        lhs.z_vertex == rhs.z_vertex;
    }
};

struct ActsSeedHash{
    std::size_t operator()(const SimSeed& s) const noexcept {
        auto getIdx = [](auto sp){
            assert(sp->sourceLinks().size() >= 1);
            auto sl = sp->sourceLinks()[0];
            return sl.template get<ActsExamples::IndexSourceLink>().index();
        };
        std::vector<std::size_t> indices;
        std::transform(s.sp().begin(), s.sp().end(), std::back_inserter(indices), getIdx);
        return indices[0] ^ indices[1] ^ indices[2];
    }
};

struct ActsSeedEquals {
    bool operator()(const SimSeed& lhs, const SimSeed& rhs) const {
        auto getIdx = [](auto sp){
            assert(sp->sourceLinks().size() >= 1);
            auto sl = sp->sourceLinks()[0];
            return sl.template get<ActsExamples::IndexSourceLink>().index();
        };
        std::vector<std::size_t> lhsIndices;
        std::transform(lhs.sp().begin(), lhs.sp().end(), std::back_inserter(lhsIndices), getIdx);
        std::vector<std::size_t> rhsIndices;
        std::transform(rhs.sp().begin(), rhs.sp().end(), std::back_inserter(rhsIndices), getIdx);
        return lhsIndices[0] == rhsIndices[0] && lhsIndices[1] == rhsIndices[1] && lhsIndices[2] == rhsIndices[2];
    }
};

/// @brief Converts a traccc seed to an acts seed.
/// @param seed the traccc seed.
/// @param spacePointConv the space point ConversionData (traccc space point -> acts space point).
/// @returns An acts seed.
template <typename T>
SimSeed convertSeed(const traccc::seed& seed, T& spacePointConv){
    return SimSeed(
        spacePointConv.indexToValue(seed.spB_link),
        spacePointConv.indexToValue(seed.spM_link),
        spacePointConv.indexToValue(seed.spT_link),
        seed.z_vertex,
        seed.weight
    );
}

/// @brief Converts a collection of traccc seeds and appends the result to the given outputContainer.
/// @param seeds the traccc seeds
/// @param SpacePointConv the spacepoint ConversionData (traccc space point -> acts space point).
/// @param outputContainer the container to put the converted space points into (an empty container is expected).
/// @returns The seed ConversionData (traccc seed -> acts seed).
template <typename T, typename allocator_t, typename output_container_t>
auto convertSeeds(std::vector<traccc::seed, allocator_t>& seeds, T& spacePointConv, output_container_t& outputContainer){
    auto fn = [&spacePointConv](auto& seed){
        return convertSeed(seed, spacePointConv);
    };
    return Util::convert<TracccSeedHash, TracccSeedEquals>(seeds, fn, outputContainer);
}

/// @brief Converts a seed to a traccc seed.
/// @param seed the seed.
/// @param spacePointConv the space point ConversionData (acts space point -> traccc space point).
/// @returns A traccc seed.
template <typename T>
traccc::seed convertSeed(const SimSeed& seed, T& spacePointConv){
    using Scalar = typename traccc::point3::value_type;
    return traccc::seed{
        spacePointConv.valueToIndex(*seed.sp()[0]),
        spacePointConv.valueToIndex(*seed.sp()[1]),
        spacePointConv.valueToIndex(*seed.sp()[2]),
        static_cast<Scalar>(seed.z()), 
        seed.seedQuality()};
}

/// @brief Converts a collection of seeds to traccc seeds and appends the result to the given outputContainer.
/// @param seeds theseeds
/// @param SpacePointConv the spacepoint ConversionData (acts space point -> traccc space point).
/// @param outputContainer the container to put the converted space points into (an empty container is expected).
/// @returns The seed ConversionData (acts seed -> traccc seed).
template <typename T, typename output_container_t>
auto convertSeeds(const SimSeedContainer& seeds, T& spacePointConv, output_container_t& outputContainer){
    auto fn = [&spacePointConv](const SimSeed& seed){
        return convertSeed(seed, spacePointConv);
    };
    return Util::convert<ActsSeedHash, ActsSeedEquals>(seeds, fn, outputContainer);
}

}  // namespace Acts::TracccPlugin
