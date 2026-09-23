// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GlobalPatternFinder.hpp"


#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"

#include "Acts/Utilities/detail/periodic.hpp"
#include "Acts/Definitions/Units.hpp"

namespace Acts::Experimental::detail {

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
template<PatternSeedSelector<Hit_t> SeedSelector_t,
         OnlyPhiHitsProvider<Hit_t, Topology_t, typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::PatternState> OnlyPhiProvider_t>
std::vector<typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::OutputPattern> 
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::findPatterns(
    const GeometryContext& gctx,
    const SearchTree_t& treeData,
    const SeedSelector_t& seedSelector,
    const OnlyPhiProvider_t& onlyPhiProvider,
    const BeamspotInfo& beamspotInfo) const 
{
    PatternStateVec patterns{findPatternsInEta(gctx, treeData, seedSelector, beamspotInfo)};
    addPhiOnlyHits(gctx, onlyPhiProvider, patterns);
    return convertToPattern(std::move(patterns));
}

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
template<PatternSeedSelector<Hit_t> SeedSelector_t>
std::vector<typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::PatternState> 
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::findPatternsInEta(
    const GeometryContext& gctx,
    const SearchTree_t& orderedSpacepoints,
    const SeedSelector_t& seedSelector,
    const BeamspotInfo& beamspotInfo) const 
{
    constexpr auto thetaIdx {Acts::toUnderlying(HitCoords::eTheta)};
    constexpr auto sectorIdx {Acts::toUnderlying(HitCoords::eSector)};

    /** Define candidate hit buffer */
    std::vector<OrderedHit> OrderedHits{};
    OrderedHits.reserve(Topology_t::nGroups * 15);

    /** Define two PatternState buffers to avoid reallocations */
    PatternStateVec startPatternBuff{}, endPatternBuff{};
    startPatternBuff.reserve(10);
    endPatternBuff.reserve(10);
    
    PatternStateVec outPatterns{};
    outPatterns.reserve(10);
    /** @brief Helper function to count existing patterns containing a hit
     *  @param hit The hit to check
     *  @param coords The coordinates of the hit
     *  @return The number of existing patterns containing the hit TO DO remove window from config */
    auto countPatterns = [&](const PatternStateVec& patterns,
                             const Hit_t& seed,
                             const SearchTree_t::coordinate_t& coords) -> unsigned {
        return std::ranges::count_if(patterns, [&](const PatternState& pattern){
            if (std::abs(pattern.patTheta - coords[thetaIdx]) > 2.*seedSelector.thetaSearchWindow(seed) ||
                !pattern.expSect.isNeighbour(
                    Sector_t{static_cast<typename Sector_t::Index_t>(coords[sectorIdx])})) {
                return false;
            }
            return pattern.isInPattern(seed);
        });
    };
    /** We try to build a pattern in eta starting from every hit in the three */
    for (const auto& [seedCoords, seedPtr] : orderedSpacepoints) {
        /** Get the seed hit */
        const Hit_t& seed {*seedPtr};
        if (!seedSelector.goodForSeeding(seed)) {
            ACTS_VERBOSE(__func__<<"() Seed hit "<<Acts::toString(*seed.spacePoint())
                <<" not good for seeding - skip.");
            continue;
        }
        ACTS_VERBOSE(__func__<<"() New seed hit "<<Acts::toString(*seed.spacePoint())
            <<" , coordinates [" << seedCoords[toUnderlying(HitCoords::eSector)] 
           << ", "<< seedCoords[toUnderlying(HitCoords::eTheta)] << "]");

        /** check how many existing patterns contain this hit */
        unsigned nExistingPatterns {countPatterns(outPatterns, seed, seedCoords)};
        if (nExistingPatterns >= m_cfg.maxSeedAttempts) {
            // Try first to resolve overlaps and re-count the number of patterns containing the seed
            outPatterns = resolveOverlaps(outPatterns);
            nExistingPatterns = countPatterns(outPatterns, seed, seedCoords);
            if (nExistingPatterns >= m_cfg.maxSeedAttempts) {
                ACTS_VERBOSE(__func__<<"() Seed has already been used in "
                    <<static_cast<int>(nExistingPatterns)
                    <<" patterns, which is above the limit - skip this seed.");
                continue;
            }   
        }
        /** Define the search range. */
        typename SearchTree_t::range_t selectRange{};
        selectRange[sectorIdx].shrink(seedCoords[sectorIdx] - 0.1, seedCoords[sectorIdx] + 0.1);
        const double thetaHalfWindow {seedSelector.thetaSearchWindow(seed)/2.};
        selectRange[thetaIdx].shrink(seedCoords[thetaIdx] - thetaHalfWindow, seedCoords[thetaIdx] + thetaHalfWindow);

        /** Search for compatible spacepoints with the seed and check if there are enough to build a pattern */
        OrderedHits.clear();
        orderedSpacepoints.rangeSearchMapDiscard(selectRange, [&](const SearchTree_t::coordinate_t& /*coords*/,
                                                                       const Hit_t* hit) {
            OrderedHits.emplace_back(hit, 0u);
        });
        if (OrderedHits.size() < m_cfg.minStripEtaLayers + m_cfg.minPrecisionLayers) {
            ACTS_VERBOSE(__func__<<"() Found "<<OrderedHits.size()<<" candidate hits, below minimum required - skip seed.");
            continue;
        }
        /** Check that the candidate hits extend at least in two layers */
        using GroupIdx = typename Topology_t::GroupIdx;
        std::array<GroupIdx, Topology_t::nGroups> groupCounts{};
        GroupIdx nValidGroups {0u};
        for (const OrderedHit& c : OrderedHits) {
            const GroupIdx groupIdx {Topology_t::groupIndex(*c)};
            if (++groupCounts[groupIdx] == m_cfg.minGroupLayers) {
                ++nValidGroups;
            }
            if (nValidGroups >= m_cfg.minGroups) {
                break;
            }
        }
        if (nValidGroups < m_cfg.minGroups) {
            ACTS_VERBOSE(__func__<<"() Found candidate hits in "<<static_cast<int>(nValidGroups)
                <<" groups, below the minimum required - skip seed.");
            continue;
        }
        /** Sort the compatible spacepoints by global logical layer */
        std::ranges::sort(OrderedHits, [&](const OrderedHit& c1, const OrderedHit& c2){
            return Topology_t::layerSorter(*c1, *c2);
        });
        /** Assign global layer number. This will avoid re-computing it many times later */
        for (typename Topology_t::LayerIdx i = 1; i < OrderedHits.size(); ++i) {
            OrderedHits[i].globLayer = OrderedHits[i - 1].globLayer + 
                !Topology_t::sameLayer(*OrderedHits[i - 1], *OrderedHits[i]);
        }
        if (OrderedHits.back().globLayer + 1u < 
                (m_cfg.minStripEtaLayers + m_cfg.minPrecisionLayers)) {
            ACTS_VERBOSE(__func__<<"() Found "<<OrderedHits.size()<<" candidate hits on "
                <<static_cast<int>(OrderedHits.back().globLayer + 1u)
                <<" layers, below the minimum required - skip this seed.");
            continue;
        }
        if (m_logger->level() <= Acts::Logging::Level::VERBOSE) {
            ACTS_VERBOSE(__func__<<"() Found "<< OrderedHits.size()<<" candidate hits: ");
            for (const auto& c : OrderedHits) {
                ACTS_VERBOSE(__func__<<"() \t**"<<c);
            }
        }

        /** Start pattern building from the seed */
        const auto seedItr {std::ranges::find_if(OrderedHits,
            [&seed](const OrderedHit& c){ return *c == seed; })};
        assert(seedItr != OrderedHits.end());
        const OrderedHit& seedCand {*seedItr};

        PatternState patternSeed{gctx, seedCand, 
            static_cast<typename Sector_t::Index_t>(seedCoords[sectorIdx]), 
            &m_cfg, m_logger.get()};
   
        /** @brief Helper function to extend a given pattern with a set of hits. We can have pattern
            *         branching when a pattern is compatible with multiple hits on the same layer,
            *         increasing the number of active patterns. For each new hit, extendPatterns will try
            *         to extend every active pattern and remove the ones not meeting continuation criteria. 
            *  @param begin The iterator pointing to the first hit to process. 
            *  @param end The iterator pointing to the end of the hit range.
            *  @param toExtend The pattern to extend.
            *  @return The vector of resulting patterns. */
        auto processHitRange = [&](const auto begin, 
                                   const auto end, 
                                   PatternState&& toExtend) -> PatternStateVec {
            startPatternBuff.clear();
            startPatternBuff.push_back(std::move(toExtend));
            
            for (auto testItr = begin; testItr != end; ++testItr) {
                const OrderedHit& testHit {*testItr};
                if (testHit.globLayer == seedCand.globLayer) {
                    continue; // skip hits on the same layer as the seed
                }
                extendPatterns(gctx, startPatternBuff, endPatternBuff, testHit, beamspotInfo);
                // Swap the buffers for the next iteration
                std::swap(startPatternBuff, endPatternBuff);
            }
            return startPatternBuff.size() > 1 
                ? resolveOverlaps(startPatternBuff) 
                : PatternStateVec{std::move(startPatternBuff.back())};
        };

        /** First search for compatible hits from the seed layer onwards */
        PatternStateVec forwardExtended {processHitRange(std::next(seedItr), OrderedHits.end(), std::move(patternSeed))};

        /** When inverting the search direction, update last inserted hit and line anchor */
        ACTS_VERBOSE(__func__<<"() Finished forward search, found "<<forwardExtended.size()<<" forward patterns, start backward search.");
        PatternStateVec backwardExtended{};
        backwardExtended.reserve(2*forwardExtended.size());
        
        for (PatternState& pat : forwardExtended) {
            ACTS_VERBOSE(__func__<<"() Start backward search for pattern "<<PatternState::detailed(pat));
            pat.moveLineAnchorHit(seedCand);
            pat.lastInsertedHit = seedCand;

            /** Then try to proceed toward the innermost layer */
            std::ranges::move(
                processHitRange(std::reverse_iterator(seedItr), OrderedHits.rend(), std::move(pat)), 
                std::back_inserter(backwardExtended)
            );
        }
        /** Ensure there are no overlaps */
        if (backwardExtended.size() > 1) {
            backwardExtended = resolveOverlaps(backwardExtended);
        }

        for (PatternState& pat : backwardExtended) {
            pat.meanNormResidual2 /= pat.nBendingLayers();
            if (!passPatternCuts(pat)) {
                continue;
            }
            ACTS_VERBOSE(__func__<<"() Add new pattern "<<PatternState::detailed(pat));
            pat.isFinalized = true;
            outPatterns.push_back(std::move(pat));
        }
    }
    ACTS_VERBOSE(__func__<<"() Found in total "<<outPatterns.size()
        <<" patterns in eta before overlap removal");
    return resolveOverlaps(outPatterns);
}
template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
void GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::extendPatterns(
    const Acts::GeometryContext& gctx,
    PatternStateVec& startPatterns,
    PatternStateVec& endPatterns,
    const OrderedHit& testHit,
    const BeamspotInfo& beamSpot) const 
{
    endPatterns.clear();
    ACTS_VERBOSE(__func__<<"() *** Test "<<testHit<<" against " 
        << startPatterns.size() << " active patterns.");

    // Compute the minimum number of missed layer hits among the active patterns, 
    // to use as reference for pruning patterns with too many missed layers. 
    auto missedLayers = [&testHit](const PatternState& pat) -> unsigned {
        return std::abs(pat.lastInsertedHit.globLayer - testHit.globLayer);
    };

    unsigned minMissedLayers {std::numeric_limits<unsigned>::max()};
    std::ranges::for_each(startPatterns, [&missedLayers, &minMissedLayers](const PatternState& pat){
        minMissedLayers = std::min(minMissedLayers, missedLayers(pat));
    });

    const bool shouldPrune {startPatterns.size() > 1 && 
        std::ranges::any_of(startPatterns, [](const PatternState& p){
            return p.nBendingLayers() > 2; })};
    
    for (auto [i, pat] : Acts::enumerate(startPatterns)) {
        if (pat.isOverlap) {
            continue;
        }
        /** Check the pattern has not already missed too many layers compared to other patterns. */
        if (Topology_t::groupIndex(*pat.lastInsertedHit) == Topology_t::groupIndex(*testHit) && 
                missedLayers(pat) > std::max(m_cfg.maxMissLayersInGroup, minMissedLayers)) {
            ACTS_VERBOSE(__func__<<"() Pattern "<<PatternState::detailed(pat)<<"\nhas missed " 
                <<static_cast<int>(missedLayers(pat))<<" layer hits, above the max allowed - abort pattern.");
            continue;
        }
        /** Prunes pattern hypotheses within groups sharing the same last-hit layer. This step reduces branching by  
         *  keeping only the best-scoring pattern within each last-hit equivalence group, while preserving all 
         *  patterns when the last-hit layer matches the reference layer (to allow further branching). */
        if (shouldPrune && pat.lastInsertedHit.globLayer != testHit.globLayer &&
            std::ranges::find_if(std::next(startPatterns.begin(), i + 1), startPatterns.end(), 
                [&] (PatternState& p){
                    if (p.lastInsertedHit != pat.lastInsertedHit || p.isOverlap) {
                        return false;
                    }
                    if (isBetter(pat, p)) {
                        ACTS_VERBOSE("extendPatterns() Pruning: "
                            <<PatternState::detailed(pat)<<"\nis BETTER than "<<PatternState::detailed(p));
                        p.isOverlap = true;
                        return false;
                    }
                    ACTS_VERBOSE("extendPatterns() Pruning: "
                        <<PatternState::detailed(p)<<"\nis BETTER than "<<PatternState::detailed(pat));
                    return true; }
            ) != startPatterns.end()) {
            continue;
        }
        /** Check angular compatibility of the test hit and the pattern */
        using LineTestDecision = PatternState::LineTestDecision;
        const auto [residual, resSigma, result] {checkLineCompatibility(gctx, pat, testHit, beamSpot)};
        switch (result) {
            case LineTestDecision::eAddHit: {
                /** TO DO: Study feasibility of loosening the criteria for low-confidence hits with OR */
                const bool lowConfidenceRes {resSigma > m_cfg.lowConfidenceResSigma && 
                                             residual / resSigma > 2.};
                if (lowConfidenceRes) {
                    ACTS_VERBOSE(__func__<<"() Low-confidence hit: residual pull "<<residual / resSigma);
                    /** If hit is compatible but with poor confidence, we create both a pattern with the hit and a pattern without the hit,
                     *  to keep also the possibility of rejecting this hit in the next iterations. First we make sure that the low-confidence
                     *  pattern is original, i.e. accumulating not seen hits */
                    if (std::ranges::any_of(endPatterns, [&testHit, &pat](const PatternState& p) {
                            return p.lastInsertedHit == testHit && 
                                   (p.prevLayerHit == pat.lastInsertedHit || 
                                        p.nBendingLayers() > (pat.nBendingLayers() + 1u)); })) {
                        ACTS_VERBOSE(__func__<<"() Forking leads to existing pattern - reject.");
                        break;
                    }
                    /** Add the new pattern to the list of next patterns */
                    endPatterns.push_back(pat);
                    endPatterns.back().addHit(gctx, testHit, residual, resSigma);
                    break;
                }
                ACTS_VERBOSE(__func__<<"() Hit compatible - add to pattern. Residual pull "<<residual / resSigma);
                pat.addHit(gctx, testHit, residual, resSigma);
                break;
            }
            case LineTestDecision::eBranchPattern: {
                /* Check first if the branched pattern already exists*/
                if (std::ranges::any_of(endPatterns, [&testHit, &pat](const PatternState& p) {
                        return p.lastInsertedHit == testHit && p.prevLayerHit == pat.prevLayerHit; })) {
                    ACTS_VERBOSE(__func__<<"() Hit compatible & on same layer of last added hit - branched pattern already exists.");
                    break;
                }
                /** Branch the pattern: we clone it and overwrite the existing hit with the test hit */
                ACTS_VERBOSE(__func__<<"() Hit compatible & on same layer of last added hit - branch pattern.");
                endPatterns.push_back(pat);
                endPatterns.back().overWriteHit(gctx, testHit, residual, resSigma);
                break;
            }
            case LineTestDecision::eRejectHit: {
                ACTS_VERBOSE(__func__<<"() Hit is not compatible with the pattern - reject hit.");
                break;
            }
        }
        endPatterns.push_back(std::move(pat));
    }
    startPatterns.clear();
};

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::LineTestRes 
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::checkLineCompatibility(
    const GeometryContext& gctx,
    PatternState& pat,
    const OrderedHit& testHit,
    const BeamspotInfo& beamSpot) const
{
    if (testHit->spacePoint()->measuresLoc0() && !pat.isPhiCompatible(gctx, *testHit)) {
        ACTS_VERBOSE(__func__<<"() Test hit phi "<<VectorHelpers::phi(testHit->globalPosition(gctx))
            <<" not compatible with "<<PatternState::brief(pat));
        return LineTestRes{};
    }
    
    /** @brief Helper function to make the result
     *  @param decision The decision for the test result if the residual is within the acceptance window 
     *  @return The test result */
    using LineTestDecision = PatternState::LineTestDecision;
    auto makeResult = [&](const LineTestDecision decision) -> LineTestRes {
        LineTestRes res {pat.computeLineResidual(gctx, testHit, beamSpot)};
        double accWindow {m_cfg.nResidualSigma * res.sigma};
        /** Loosen the window when we use the beamspot or when we are looking for hits in a new group, as
            *  the straight line approximation becomes less accurate on large distances. TO DO: investigate this further */
        if (pat.useBeamspot || 
            Topology_t::groupIndex(*testHit) != Topology_t::groupIndex(*pat.lastInsertedHit) ||
            (Topology_t::groupIndex(*testHit) != Topology_t::groupIndex(*pat.prevLayerHit) 
                && testHit.globLayer == pat.lastInsertedHit.globLayer)) {
            accWindow *= 2.;
        }
        if (res.residual < accWindow) {
            res.result = decision;
        }
        return res;
    };

    if(testHit.globLayer != pat.lastInsertedHit.globLayer) {
        pat.updateLineParameters(gctx, beamSpot);
        return makeResult(LineTestDecision::eAddHit);
    }
    if (testHit == pat.lastInsertedHit) {
        ACTS_VERBOSE(__func__<<"() Test hit is the same as last inserted hit - reject.");
        return LineTestRes{};
    }
    if (pat.lineAnchorHit.globLayer == pat.lastInsertedHit.globLayer) {
        ACTS_VERBOSE(__func__<<"() Test hit on same layer as seed with no prior hits - reject.");
        return LineTestRes{};
    }
    return makeResult(LineTestDecision::eBranchPattern);
}

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
bool GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::passPatternCuts(
    const PatternState& pat) const 
{
    /** Check that the pattern meets the minimum requirements for trigger and precision layers */
    if (pat.nTriggerLayers < m_cfg.minStripEtaLayers || 
        pat.nPrecisionLayers < m_cfg.minPrecisionLayers ||
        std::ranges::count_if(pat.hitsPerGroup, 
            [this](const auto& hits) { return hits.size() >= m_cfg.minGroupLayers; }) < 2) {
        ACTS_VERBOSE(__func__<<"() Pattern " << PatternState::detailed(pat)
            << "\ndoes not meet minimum layer requirements - reject.");
        return false;
    }
    /** Check requirement on the residual */
    if (pat.meanNormResidual2 > m_cfg.meanNormRes2Cut) {
        ACTS_VERBOSE(__func__<<"() Pattern " << PatternState::detailed(pat)
            << "\ndoes not meet the mean norm residual2 cut - reject.");
        return false;
    }
    return true;
}

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
bool 
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::isBetter(
    const PatternState& a, 
    const PatternState& b) 
{
    const double resA {a.getMeanResidual2()};
    const double resB {b.getMeanResidual2()};
    const double resDiff {
        std::abs(resA - resB) / std::max(resA, resB)
    };
    const int nLayerDiff {a.nBendingLayers() - b.nBendingLayers()};
    const int nPrecLayDiff {a.nPrecisionLayers - b.nPrecisionLayers};

    /** For patterns that differ by 1–2 layers, don't sacrifice fit quality  
     *  unless the extra layers are genuinely comparable. For a ≥3-layer difference,  
     *  the multiplicity advantage is strong enough to dominate. */
    if ((nLayerDiff == 0 && nPrecLayDiff == 0) ||
        (std::abs(nLayerDiff) < 3 && resDiff > 0.1)) {
        return resA < resB;
    }
    if (nLayerDiff == 0) {
        return nPrecLayDiff > 0;
    }
    return nLayerDiff > 0;
}

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::PatternStateVec
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::resolveOverlaps(
    PatternStateVec& toResolve) const 
{
    ACTS_VERBOSE(__func__<<"() Resolving overlaps among "
        <<toResolve.size()<<" patterns.");
    PatternStateVec outputPatterns{};
    outputPatterns.reserve(toResolve.size());
    /** Check if two patterns overlap in space */
    auto areOverlapping = [this](const PatternState& a, const PatternState& b) {
        /** Check first the geometrical overlap */
        if(!a.expSect.isNeighbour(b.expSect)) {
            return false;
        }
        /** Check the angular difference between the seed hits */
        if (std::abs(a.patTheta - b.patTheta) > 2.*m_cfg.maxThetaSizeOverlap) {
            return false;
        }
        if (a.nPhiLayers > 0 && b.nPhiLayers > 0) {
            using namespace Acts::UnitLiterals;
            if (std::abs(Acts::detail::difference_periodic(
                    a.patPhi, b.patPhi, 2. * std::numbers::pi)) 
                        > 5._degree) {
                return false;
            }
        } else if (a.nPhiLayers > 0) {
            if (!b.expSect.insideSector(a.patPhi)) {
                return false;
            }
        } else if (b.nPhiLayers > 0) {
            if (!a.expSect.insideSector(b.patPhi)) {
                return false;
            }
        }
        /** If we reach here, the patterns can overlap geometrically, so check the hit content */
        std::size_t nSharedHits{0}, nSharedGroups{0};
        for (std::size_t st{0u}; st < Topology_t::nGroups; ++st) {
            const auto& hitsA {a.hitsPerGroup[st]};
            const auto& hitsB {b.hitsPerGroup[st]};
            if (hitsA.empty() || hitsB.empty()) {
                continue;
            }
            const std::size_t nSharedInGroup = std::ranges::count_if(hitsA, [&](const OrderedHit& hitA){
                return std::ranges::any_of(hitsB, [&hitA](const OrderedHit& hitB) {
                    return hitA->spacePoint()->primaryMeasurement() == 
                           hitB->spacePoint()->primaryMeasurement();
                });
            });
            nSharedHits += nSharedInGroup;
            if (nSharedInGroup >= m_cfg.minGroupLayers) {
                nSharedGroups++;
            }
        }
        /** Overlap if more than 50% of the hits of the smaller pattern are shared */
        const std::size_t minHits {std::min(a.nBendingLayers(), 
                                            b.nBendingLayers())};
        const std::size_t minGroups {std::min(a.countGroups(/*onlyGoodGroups=*/ true), 
                                              b.countGroups(/*onlyGoodGroups=*/ true))};
        return nSharedHits >= 0.5 *minHits && 
               nSharedGroups >= std::min(2ul, minGroups);
    };
    /** Determine best pattern */
    auto isBetterOverlap = [&](const PatternState& a, const PatternState& b) {
        const int nGoodGroupDiff {a.countGroups(/*onlyGoodGroups=*/ true) - 
                                  b.countGroups(/*onlyGoodGroups=*/ true)};
        if (nGoodGroupDiff != 0) {
            return nGoodGroupDiff > 0;
        }
        return isBetter(a,b);
    };

    for (auto it = toResolve.begin(); it != toResolve.end(); ++it) {
        if (it->isOverlap) {
            continue;
        }
        for (auto jt = std::next(it); jt != toResolve.end(); ++jt) {
            if (jt->isOverlap || !areOverlapping(*it, *jt)) {
                continue;
            }
            if (isBetterOverlap(*it, *jt)) {
                ACTS_VERBOSE(__func__<<"() Pattern "<<PatternState::detailed(*it)
                    <<"\nis BETTER than "<<PatternState::detailed(*jt));
                jt->isOverlap = true;
            } else {
                it->isOverlap = true;
                ACTS_VERBOSE(__func__<<"() Pattern "<<PatternState::detailed(*jt)
                    <<"\nis BETTER than "<<PatternState::detailed(*it));
                break;
            }
        }
        if (!it->isOverlap) {
            outputPatterns.push_back( std::move(*it));
        }
    }
    ACTS_VERBOSE(__func__<<"() Patterns surviving overlap removal: "<< outputPatterns.size());
    return outputPatterns;
}

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
template<OnlyPhiHitsProvider<Hit_t, Topology_t, 
                             typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::PatternState> OnlyPhiProvider_t>
void GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::addPhiOnlyHits(
    const GeometryContext& gctx,
    const OnlyPhiProvider_t& onlyPhiProvider,
    PatternStateVec& patterns) const 
{
    auto computePatternLineInGroup = [&](PatternState& pat, 
                                         const Hit_t& phiHit) -> bool {
        const auto phiHitGroup {Topology_t::groupIndex(phiHit)};
        const auto& groupHits {pat.hitsPerGroup[phiHitGroup]};

        /** We use useBeamspot as a flag to indicate whether the pattern line has been determined successfully */
        pat.useBeamspot = true;

        if (groupHits.size() > 1) {
            // if we have >= 2 eta hits in the group, we use the furthestmost to define the pattern line
            const auto [minIt, maxIt] {std::ranges::minmax_element(groupHits, {},
                [](const OrderedHit& c){ return c.globLayer; })};
            pat.lineAnchorHit = *minIt;
            pat.lastInsertedHit = *maxIt;
            pat.updateLineParameters(gctx, BeamspotInfo{});
        }
        if (pat.useBeamspot && !groupHits.empty()) {
            // if we have only one eta hit or the layer separation is too small, to find the second hit 
            // we use the functionality of anchor hit
            pat.moveLineAnchorHit(groupHits.front());
            pat.lastInsertedHit = *std::ranges::max_element(groupHits, {},
                [&](const OrderedHit& c){ 
                    return (pat.projToPhiPlane(gctx, *c) - 
                            pat.projToPhiPlane(gctx, *pat.lineAnchorHit)).norm(); }); 
            pat.updateLineParameters(gctx, BeamspotInfo{});
        }
        if (pat.useBeamspot) {
            /** Define the closest hits in the upper and lewer groups */
            const OrderedHit *lowerHit{nullptr}, *upperHit{nullptr};
            for (const auto& [group, hits] : Acts::enumerate(pat.hitsPerGroup)) {
                if (hits.empty() || group == phiHitGroup) {
                    continue;
                }
                if (Topology_t::layerSorter(*hits.front(), phiHit)) {
                    const OrderedHit& lowerInGroup {*std::ranges::max_element(hits, {},
                        [](const OrderedHit& c){ return c.globLayer; })};
                    if (!lowerHit || lowerInGroup.globLayer > lowerHit->globLayer) {
                        lowerHit = &lowerInGroup;
                    }
                } else {
                    const OrderedHit& upperInGroup {*std::ranges::min_element(hits, {},
                        [](const OrderedHit& c){ return c.globLayer; })};
                    if (!upperHit || upperInGroup.globLayer < upperHit->globLayer) {
                        upperHit = &upperInGroup;
                    }
                }
            }
            if (!lowerHit || !upperHit) {
                return false;
            }
            pat.lineAnchorHit = *lowerHit;
            pat.lastInsertedHit = *upperHit;
            pat.updateLineParameters(gctx, BeamspotInfo{});
        }
        return !pat.useBeamspot;
    };

    PatternStateVec survivingPatterns{};
    survivingPatterns.reserve(patterns.size());
    for (PatternState& pat : patterns) {
        /** We look for phi-only hits in the buckets associated with the pattern */
        ACTS_VERBOSE(__func__<<"() Search for phi-only hits for pattern: " << PatternState::brief(pat));

        std::optional<typename Topology_t::GroupIdx> patterLineGroup{std::nullopt};

        auto projOntoPhiPlane = [&pat](const Vector3& pos) -> Vector3 {
            return pos - pos.dot(pat.bendPlaneNorm) * pat.bendPlaneNorm;
        };
        auto phiPull = [&](const Hit_t& hit) {
            return std::abs(Acts::detail::difference_periodic(
                        VectorHelpers::phi(hit.globalPosition(gctx)), pat.patPhi, 2. * std::numbers::pi)
                    ) / std::sqrt(hit.phiVariance(gctx));
        };

        std::array<std::vector<Hit_t>, Topology_t::nGroups> phihitsPerGroup{
            onlyPhiProvider.getPhiOnlyHits(pat, gctx)};

        for (auto [group, bucket] : Acts::enumerate(phihitsPerGroup)) {
            
            for (Hit_t& newHit : bucket) {
                ACTS_VERBOSE(__func__<<"() *** Test phi-only hit "
                    <<Acts::toString(*newHit.spacePoint())<<" in group "<<static_cast<int>(group));
                
                /** Reject hits from a layer that already contains a phi hit */
                const auto& groupHits{pat.hitsPerGroup[group]};
                if (std::ranges::any_of(groupHits, 
                        [&](const OrderedHit& h){
                            return h->spacePoint()->measuresLoc0() && 
                                Topology_t::sameLayer(*h, newHit); })) {
                    ACTS_VERBOSE(__func__<<"() Found eta hit measuring also phi on same layer - skip hit.");
                    continue;
                } 

                if (!pat.isPhiCompatible(gctx, newHit)) {
                    ACTS_VERBOSE(__func__<<"() Phi-only hit not compatible");
                    continue;
                }

                if (!patterLineGroup || *patterLineGroup != group) {
                    if (!computePatternLineInGroup(pat, newHit)) {
                        ACTS_VERBOSE(__func__<<"() Invalid projection model for group "<<group<<" - skip hit.");
                        continue;
                    }
                    patterLineGroup = group;
                }
                const Vector3 newHitPos {newHit.globalPosition(gctx)};
                using SeedingAux = Acts::Experimental::detail::CompSpacePointAuxiliaries;
                const double stripHalfLength {
                    std::sqrt(newHit.spacePoint()->covariance()[Acts::toUnderlying(SeedingAux::ResidualIdx::bending)])};
                const Vector3 sensorDir {newHit.globalSensorDirection(gctx)};
                const Vector3 stripLow {
                    projOntoPhiPlane(newHitPos - stripHalfLength * sensorDir)};
                const Vector3 stripHigh {
                    projOntoPhiPlane(newHitPos + stripHalfLength * sensorDir)};
                const Vector3 stripDirOnPlane {(stripHigh - stripLow).normalized()};

                const double stripIntersect {Acts::detail::LineHelper::lineIntersect<3>(
                    pat.linePos, pat.lineDir, stripLow, stripDirOnPlane).pathLength()};
                const double stripProjLength {(stripHigh - stripLow).norm()};

                ACTS_VERBOSE(__func__<<"() Intersect distance from lower strip edge: "
                    <<stripIntersect<<", proj strip length: "<<stripProjLength);
                
                using namespace Acts::UnitLiterals;
                constexpr double margin {10._mm};
                if (stripIntersect < -margin || stripIntersect > (stripProjLength + margin)) {
                    ACTS_VERBOSE(__func__<<"() The pattern falls outside the test hit strip in eta - skip hit.");
                    continue;
                }

                /** Reject hits from a layer that already contains a phi hit */
                if (const auto it = std::ranges::find_if(pat.phiOnlyHits, 
                        [&](const Hit_t& h){
                            return Topology_t::sameLayer(h, newHit); });
                    it != pat.phiOnlyHits.end()) {
                    if (phiPull(*it) < phiPull(newHit)) {
                        ACTS_VERBOSE(__func__<<"() A phi hit with better pull exists on same layer - skip hit.");
                        continue;
                    }
                    *it = std::move(newHit);
                } else {
                    pat.phiOnlyHits.push_back(std::move(newHit));
                    pat.nPhiLayers++;
                }
                pat.updatePatternPhi(gctx);
            }
        }
        if (pat.nPhiLayers < m_cfg.minPhiLayers) {
            ACTS_VERBOSE(__func__<<"() Pattern "<<PatternState::detailed(pat)
                <<" has only "<<static_cast<int>(pat.nPhiLayers)
                <<" phi layers, below the minimum required - reject this pattern.");
            continue;
        }
        survivingPatterns.push_back(std::move(pat));
    }
    std::swap(patterns, survivingPatterns);
}


template<GlobPatFinderHit Hit_t,
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

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::OutputPattern 
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::convertToPattern(
    PatternState&& candidate) const 
{
    OutputPattern output{static_cast<typename Sector_t::Index_t>(candidate.expSect.sector())};
    for (typename Topology_t::GroupIdx g = 0u; g < Topology_t::nGroups; ++g) {
        output.hitsPerGroup[g].reserve(candidate.hitsPerGroup[g].size());
        std::ranges::transform(candidate.hitsPerGroup[g], std::back_inserter(output.hitsPerGroup[g]), 
            [](const OrderedHit& hit) { return hit.hitPtr; });
    }
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

template<GlobPatFinderHit Hit_t,
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