// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GlobalPatternFinder.hpp"

namespace Acts::Experimental::detail {

template<HitPayload Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::GlobalPatternFinder(
    Config&& config,
    std::unique_ptr<const Logger> logger) :
    m_cfg{std::move(config)},
    m_logger{std::move(logger)} 
{
    static_assert(std::is_move_assignable_v<PatternState>);
    static_assert(std::is_move_constructible_v<PatternState>);
    static_assert(std::is_copy_assignable_v<PatternState>);
    static_assert(std::is_copy_constructible_v<PatternState>);
    static_assert(std::is_nothrow_move_constructible_v<PatternState>);
    static_assert(std::is_nothrow_move_assignable_v<PatternState>);
}

template<HitPayload Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::OutputPattern 
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::convertToPattern(
    PatternState&& candidate) const 
{
    OutputPattern output{};
    output.hits = std::move(candidate.hitsPerGroup);
    output.phiOnlyHits = std::move(candidate.phiOnlyHits);
    output.meanNormResidual2 = candidate.meanNormResidual2;
    output.patTheta = candidate.patTheta;
    output.patPhi = candidate.patPhi;
    output.expSect = Sector_t{candidate.expSect.sector()};
    output.nPrecisionLayers = candidate.nPrecisionLayers;
    output.nTriggerLayers = candidate.nTriggerLayers;
    output.nPhiLayers = candidate.nPhiLayers;

    return output;
}

template<HitPayload Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
std::vector<typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::OutputPattern>
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::convertToPattern(
    PatternStateVec&& candidates) const 
{
    std::vector<OutputPattern> outPatterns{};
    outPatterns.reserve(candidates.size());

    for (PatternState& pat : candidates) {
        outPatterns.push_back(convertToPattern(std::move(pat)));
    }
    return outPatterns;
}
}