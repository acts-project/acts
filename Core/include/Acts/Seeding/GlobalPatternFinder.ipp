// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GlobalPatternFinder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

namespace Acts::Experimental::detail {

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
template <PatternSeedSelector<Hit_t, Topology_t> SeedSelector_t,
          OnlyPhiHitsProvider<Hit_t, Topology_t, Sector_t> OnlyPhiProvider_t>
std::vector<
    typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::OutputPattern>
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::findPatterns(
    const GeometryContext& gctx, const SearchTree_t& treeData,
    const SeedSelector_t& seedSelector,
    const OnlyPhiProvider_t& onlyPhiProvider,
    const BeamspotInfo& beamspotInfo) const {
  PatternStateAuxVec patterns{
      findPatternsInEta(gctx, treeData, seedSelector, beamspotInfo)};
  if (m_cfg.searchOnlyPhiHits) {
    addPhiOnlyHits(gctx, beamspotInfo, onlyPhiProvider, patterns);
  }
  return convertToPattern(std::move(patterns));
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
template <PatternSeedSelector<Hit_t, Topology_t> SeedSelector_t>
std::vector<
    typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::PatternStateAux>
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::findPatternsInEta(
    const GeometryContext& gctx, const SearchTree_t& orderedSpacepoints,
    const SeedSelector_t& seedSelector,
    const BeamspotInfo& beamspotInfo) const {
  constexpr auto thetaIdx{Acts::toUnderlying(HitCoords::eTheta)};
  constexpr auto sectorIdx{Acts::toUnderlying(HitCoords::eSector)};

  /** Define candidate hit buffer */
  std::vector<OrderedHit> orderedHits{};
  orderedHits.reserve(Topology_t::nGroups * 15);

  /** Define two PatternStateAux buffers to avoid reallocations */
  PatternStateAuxVec startPatternBuff{}, endPatternBuff{};
  startPatternBuff.reserve(10);
  endPatternBuff.reserve(10);

  PatternStateAuxVec outPatterns{};
  outPatterns.reserve(10);

  /** We try to build a pattern in eta starting from every hit in the tree */
  for (const auto& [seedCoords, seedPtr] : orderedSpacepoints) {
    /** Get the seed hit */
    const Hit_t& seed{*seedPtr};
    if (!seedSelector.goodForSeeding(seed)) {
      continue;
    }
    ACTS_VERBOSE(__func__ << "() New seed hit "
                          << Acts::toString(*seed.spacePoint())
                          << " , coordinates [" << seedCoords[sectorIdx] << ", "
                          << seedCoords[thetaIdx] << "]");

    /** Check how many existing patterns contain this hit */
    if (countPatterns(outPatterns, seed, seedCoords) >= m_cfg.maxSeedAttempts) {
      ACTS_VERBOSE(__func__
                   << "() Seed has already been used in "
                   << static_cast<int>(
                          countPatterns(outPatterns, seed, seedCoords))
                   << " patterns, which is above the limit - skip this seed.");
      continue;
    }
    /** Define the search range. */
    const auto& windowPerGroup{seedSelector.thetaSearchWindow(seed)};
    const double thetaHalfWindow{std::ranges::max(windowPerGroup)};

    typename SearchTree_t::range_t selectRange{};
    // The sector coordinate is an integer index: +-0.1 selects only the seed
    // sector
    selectRange[sectorIdx].shrink(seedCoords[sectorIdx] - 0.1,
                                  seedCoords[sectorIdx] + 0.1);
    selectRange[thetaIdx].shrink(seedCoords[thetaIdx] - thetaHalfWindow,
                                 seedCoords[thetaIdx] + thetaHalfWindow);

    /** Search for compatible spacepoints with the seed and check if there are
     * enough to build a pattern */
    orderedHits.clear();
    orderedSpacepoints.rangeSearchMapDiscard(
        selectRange,
        [&](const SearchTree_t::coordinate_t& coords, const Hit_t* hit) {
          if (std::abs(coords[thetaIdx] - seedCoords[thetaIdx]) >
              windowPerGroup[Topology_t::groupIndex(*hit)]) {
            return;  // skip hits outside the theta search window
          }
          orderedHits.emplace_back(hit, 0u);
        });
    if (orderedHits.size() <
        m_cfg.minStripEtaLayers + m_cfg.minPrecisionLayers) {
      ACTS_VERBOSE(__func__
                   << "() Found " << orderedHits.size()
                   << " candidate hits, below minimum required - skip seed.");
      continue;
    }
    /** Check that at least minGroups groups contain minGroupLayers candidate
     * hits each */
    using GroupIdx = typename Topology_t::GroupIdx;
    std::array<GroupIdx, Topology_t::nGroups> groupCounts{};
    GroupIdx nValidGroups{0u};
    for (const OrderedHit& c : orderedHits) {
      const GroupIdx groupIdx{Topology_t::groupIndex(*c)};
      if (++groupCounts[groupIdx] == m_cfg.minGroupLayers) {
        ++nValidGroups;
      }
      if (nValidGroups >= m_cfg.minGroups) {
        break;
      }
    }
    if (nValidGroups < m_cfg.minGroups) {
      ACTS_VERBOSE(__func__
                   << "() Found candidate hits in "
                   << static_cast<int>(nValidGroups)
                   << " groups, below the minimum required - skip seed.");
      continue;
    }
    /** Sort the compatible spacepoints by global logical layer */
    std::ranges::sort(orderedHits,
                      [&](const OrderedHit& c1, const OrderedHit& c2) {
                        return Topology_t::layerSorter(*c1, *c2);
                      });
    /** Assign global layer number. This will avoid re-computing it many times
     * later */
    for (std::size_t i = 1; i < orderedHits.size(); ++i) {
      orderedHits[i].globLayer =
          orderedHits[i - 1].globLayer +
          !Topology_t::sameLayer(*orderedHits[i - 1], *orderedHits[i]);
    }
    if (orderedHits.back().globLayer + 1u <
        (m_cfg.minStripEtaLayers + m_cfg.minPrecisionLayers)) {
      ACTS_VERBOSE(__func__
                   << "() Found " << orderedHits.size() << " candidate hits on "
                   << static_cast<int>(orderedHits.back().globLayer + 1u)
                   << " layers, below the minimum required - skip this seed.");
      continue;
    }
    if (m_logger->level() <= Acts::Logging::Level::VERBOSE) {
      ACTS_VERBOSE(__func__ << "() Found " << orderedHits.size()
                            << " candidate hits: ");
      for (const auto& c : orderedHits) {
        ACTS_VERBOSE(__func__ << "() \t**" << c);
      }
    }

    /** Start pattern building from the seed */
    const auto seedItr{std::ranges::find_if(
        orderedHits,
        [&seedPtr](const OrderedHit& c) { return c.hitPtr == seedPtr; })};
    assert(seedItr != orderedHits.end());
    const OrderedHit& seedCand{*seedItr};

    PatternStateAux patternSeed{
        gctx, seedCand,
        static_cast<typename Sector_t::Index_t>(seedCoords[sectorIdx]), &m_cfg,
        m_logger.get()};

    /** @brief Helper function to extend a given pattern with a set of hits. We can have pattern
     *         branching when a pattern is compatible with multiple hits on the
     * same layer, increasing the number of active patterns. For each new hit,
     * extendPatterns will try to extend every active pattern and remove the
     * ones not meeting continuation criteria.
     *  @param begin The iterator pointing to the first hit to process.
     *  @param end The iterator pointing to the end of the hit range.
     *  @param isForward Whether the search is forward or backward in the hit range.
     *  @param toExtend The pattern to extend.
     *  @return The vector of resulting patterns. */
    auto processHitRange =
        [&](const auto begin, const auto end, const bool isForward,
            PatternStateAux&& toExtend) -> PatternStateAuxVec {
      startPatternBuff.clear();
      startPatternBuff.push_back(std::move(toExtend));

      for (auto testItr = begin; testItr != end; ++testItr) {
        const OrderedHit& testHit{*testItr};
        if (testHit.globLayer == seedCand.globLayer) {
          continue;  // skip hits on the same layer as the seed
        }
        extendPatterns(gctx, startPatternBuff, endPatternBuff, testHit,
                       beamspotInfo);
        // Swap the buffers for the next iteration
        std::swap(startPatternBuff, endPatternBuff);
      }
      PatternStateAuxVec result{};
      result.reserve(startPatternBuff.size());

      for (PatternStateAux& pat : startPatternBuff) {
        if (pat.isOverlap) {
          continue;
        }
        if (isForward) {
          // We need to update the last inserted hit and line anchor when the
          // search direction is inverted
          pat.lastInsertedHit = seedCand;
          pat.moveLineAnchorHit(seedCand);
        } else if (!passPatternCuts(pat)) {
          continue;
        }
        pat.syncState(gctx, beamspotInfo);
        // Resolve patterns produced from the same pattern seed
        insertIfBest(gctx, result, std::move(pat));
      }
      return result;
    };

    /** First search for compatible hits from the seed layer onwards */
    PatternStateAuxVec forwardExtended{
        processHitRange(std::next(seedItr), orderedHits.end(),
                        /*isForward*/ true, std::move(patternSeed))};

    /** Then search backward from the seed layer, starting from each forward
     * pattern */
    ACTS_VERBOSE(__func__ << "() Finished forward search, found "
                          << forwardExtended.size()
                          << " forward patterns, start backward search.");
    PatternStateAuxVec backwardExtended{};
    backwardExtended.reserve(2 * forwardExtended.size());

    for (PatternStateAux& pat : forwardExtended) {
      ACTS_VERBOSE(__func__ << "() Start backward search for pattern "
                            << PatternStateAux::detailed(pat));

      for (PatternStateAux& backwardCand :
           processHitRange(std::reverse_iterator(seedItr), orderedHits.rend(),
                           /*isForward*/ false, std::move(pat))) {
        // Resolve siblings produced by different backward searches (one for
        // each forward pattern)
        insertIfBest(gctx, backwardExtended, std::move(backwardCand));
      }
    }
    // Siblings are now mutually non-overlapping: insert the survivors into the
    // output
    for (PatternStateAux& pat : backwardExtended) {
      if (insertIfBest(gctx, outPatterns, std::move(pat))) {
        ACTS_VERBOSE(__func__ << "() Added new pattern "
                              << PatternStateAux::detailed(outPatterns.back()));
      } else {
        ACTS_VERBOSE(
            __func__
            << "() Pattern " << PatternStateAux::detailed(pat)
            << "\nnot added to output patterns because better patterns exist.");
      }
    }
  }
  ACTS_VERBOSE(__func__ << "() Found in total " << outPatterns.size()
                        << " patterns in eta.");
  return outPatterns;
}
template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
void GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::extendPatterns(
    const Acts::GeometryContext& gctx, PatternStateAuxVec& startPatterns,
    PatternStateAuxVec& endPatterns, const OrderedHit& testHit,
    const BeamspotInfo& beamSpot) const {
  endPatterns.clear();
  ACTS_VERBOSE(__func__ << "() *** Test " << testHit << " against "
                        << startPatterns.size() << " active patterns.");

  // Compute the minimum number of missed layer hits among the active patterns,
  // to use as reference for pruning patterns with too many missed layers.
  auto missedLayers = [&testHit](const PatternStateAux& pat) -> unsigned {
    return layerDistance(pat.lastInsertedHit.globLayer, testHit.globLayer);
  };

  unsigned minMissedLayers{std::numeric_limits<unsigned>::max()};
  std::ranges::for_each(startPatterns, [&missedLayers, &minMissedLayers](
                                           const PatternStateAux& pat) {
    minMissedLayers = std::min(minMissedLayers, missedLayers(pat));
  });

  const bool shouldPrune{
      startPatterns.size() > 1 &&
      std::ranges::any_of(startPatterns, [](const PatternStateAux& p) {
        return p.nBendingLayers() > 2;
      })};

  for (auto [i, pat] : Acts::enumerate(startPatterns)) {
    if (pat.isOverlap) {
      continue;
    }
    /** Check the pattern has not already missed too many layers compared to
     * other patterns. */
    if (Topology_t::groupIndex(*pat.lastInsertedHit) ==
            Topology_t::groupIndex(*testHit) &&
        missedLayers(pat) >
            std::max(m_cfg.maxMissLayersInGroup, minMissedLayers)) {
      ACTS_VERBOSE(__func__
                   << "() Pattern " << PatternStateAux::detailed(pat)
                   << "\nhas missed " << static_cast<int>(missedLayers(pat))
                   << " layer hits, above the max allowed - abort pattern.");
      continue;
    }
    /** Prunes pattern hypotheses within groups sharing the same last-hit layer.
     * This step reduces branching by keeping only the best-scoring pattern
     * within each last-hit equivalence group, while preserving all patterns
     * when the last-hit layer matches the reference layer (to allow further
     * branching). */
    if (shouldPrune && pat.lastInsertedHit.globLayer != testHit.globLayer &&
        std::ranges::find_if(std::next(startPatterns.begin(), i + 1),
                             startPatterns.end(), [&](PatternStateAux& p) {
                               if (p.lastInsertedHit != pat.lastInsertedHit ||
                                   p.isOverlap) {
                                 return false;
                               }
                               if (isBetter(pat, p)) {
                                 ACTS_VERBOSE("extendPatterns() Pruning: "
                                              << PatternStateAux::detailed(pat)
                                              << "\nis BETTER than "
                                              << PatternStateAux::detailed(p));
                                 p.isOverlap = true;
                                 return false;
                               }
                               ACTS_VERBOSE("extendPatterns() Pruning: "
                                            << PatternStateAux::detailed(p)
                                            << "\nis BETTER than "
                                            << PatternStateAux::detailed(pat));
                               return true;
                             }) != startPatterns.end()) {
      continue;
    }
    /** Check angular compatibility of the test hit and the pattern */
    ACTS_VERBOSE(__func__ << "() " << PatternStateAux::brief(pat));
    const auto [decision,
                testRes]{checkLineCompatibility(gctx, pat, testHit, beamSpot)};
    switch (decision) {
      case LineTestDecision::eAddHit: {
        /** TO DO: Study feasibility of loosening the criteria for
         * low-confidence hits with OR */
        const bool lowConfidenceRes{testRes.sigma >
                                        m_cfg.lowConfidenceResSigma &&
                                    testRes.residual / testRes.sigma > 2.};
        if (lowConfidenceRes) {
          /** If hit is compatible but with poor confidence, we create both a
           * pattern with the hit and a pattern without the hit, to keep also
           * the possibility of rejecting this hit in the next iterations. First
           * we make sure that the low-confidence pattern is original, i.e.
           * accumulating not seen hits */
          if (std::ranges::any_of(
                  endPatterns, [&testHit, &pat](const PatternStateAux& p) {
                    return p.lastInsertedHit == testHit &&
                           (p.prevLayerHit == pat.lastInsertedHit ||
                            p.nBendingLayers() > (pat.nBendingLayers() + 1u));
                  })) {
            ACTS_VERBOSE(__func__ << "() Low-confidence hit: Forking leads to "
                                     "existing pattern - reject.");
            break;
          }
          /** Add the new pattern to the list of next patterns */
          ACTS_VERBOSE(__func__
                       << "() Low-confidence hit - forking the pattern.");
          endPatterns.push_back(pat);
          endPatterns.back().addHit(testHit, testRes);
          break;
        }
        ACTS_VERBOSE(__func__ << "() Hit compatible - add to pattern.");
        pat.addHit(testHit, testRes);
        break;
      }
      case LineTestDecision::eBranchPattern: {
        /** Check first if the branched pattern already exists */
        if (std::ranges::any_of(endPatterns,
                                [&testHit, &pat](const PatternStateAux& p) {
                                  return p.lastInsertedHit == testHit &&
                                         p.prevLayerHit == pat.prevLayerHit;
                                })) {
          ACTS_VERBOSE(__func__ << "() Hit compatible & on same layer - "
                                   "branched pattern already exists.");
          break;
        }
        /** Branch the pattern: we clone it and overwrite the existing hit with
         * the test hit */
        ACTS_VERBOSE(__func__
                     << "() Hit compatible & on same layer - branch pattern.");
        endPatterns.push_back(pat);
        endPatterns.back().overWriteHit(testHit, testRes);
        break;
      }
      case LineTestDecision::eReplaceHit: {
        ACTS_VERBOSE(
            __func__
            << "() Hit compatible & on same layer - replace existing hit.");
        pat.overWriteHit(testHit, testRes);
        break;
      }
      case LineTestDecision::eRejectHit: {
        ACTS_VERBOSE(__func__ << "() Hit is not compatible - reject hit.");
        break;
      }
    }
    endPatterns.push_back(std::move(pat));
  }
  startPatterns.clear();
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::LineTestResult_t
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::checkLineCompatibility(
    const GeometryContext& gctx, PatternStateAux& pat,
    const OrderedHit& testHit, const BeamspotInfo& beamSpot) const {
  const bool searchInNewLayer{testHit.globLayer !=
                              pat.lastInsertedHit.globLayer};
  if (searchInNewLayer) {
    pat.syncState(gctx, beamSpot);
  }
  if (testHit->spacePoint()->measuresLoc0() &&
      !isPhiCompatible(gctx, pat, *testHit)) {
    return std::make_pair(LineTestDecision::eRejectHit, LineTestRes{});
  }

  /** @brief Helper function to make the result
   *  @param decision The decision for the test result if the residual is within the acceptance window
   *  @return The test result */
  auto makeResult = [&](const LineTestDecision decision) -> LineTestResult_t {
    LineTestRes res{pat.computeLineResidual(gctx, testHit, beamSpot)};
    double accWindow{(res.nDof == 1u ? m_cfg.nResidualSigma : m_cut2Dof) *
                     res.sigma};
    /** Loosen the window when we use the beamspot or when we are looking for
     * hits in a new group, as the straight line approximation becomes less
     * accurate on large distances. TO DO: investigate this further */
    const auto testHitGroup{Topology_t::groupIndex(*testHit)};
    const bool searchInNewGroup{testHitGroup !=
                                Topology_t::groupIndex(*pat.lastInsertedHit)};
    if (pat.useBeamspot || searchInNewGroup ||
        (testHitGroup != Topology_t::groupIndex(*pat.prevLayerHit) &&
         !searchInNewLayer)) {
      accWindow *= m_cfg.loosenAcceptanceFactor;
    }
    if (res.residual < accWindow) {
      return std::make_pair(decision, res);
    }
    return std::make_pair(LineTestDecision::eRejectHit, res);
  };

  if (searchInNewLayer) {
    return makeResult(LineTestDecision::eAddHit);
  }
  if (testHit == pat.lastInsertedHit) {
    ACTS_VERBOSE(
        __func__ << "() Test hit is the same as last inserted hit - reject.");
    return std::make_pair(LineTestDecision::eRejectHit, LineTestRes{});
  }
  if (pat.nBendingLayers() < 2u) {
    ACTS_VERBOSE(
        __func__
        << "() Test hit on same layer as seed with no prior hits - reject.");
    return std::make_pair(LineTestDecision::eRejectHit, LineTestRes{});
  }
  /** Same layer, one precision and one non-precision hit: decide whether they
   * are the same crossing (keep only the precision hit) or whether both
   * hypotheses must be kept (branch). */
  const bool testPrec{testHit->isPrecision()};
  const bool lastPrec{pat.lastInsertedHit->isPrecision()};
  if (testPrec != lastPrec) {
    ACTS_VERBOSE(__func__ << "() Test hit on same layer as last inserted hit "
                             "with different precision - check compatibility.");
    LineTestResult_t result{makeResult(LineTestDecision::eBranchPattern)};
    // If the test hit is not compatible, reject it
    if (result.first == LineTestDecision::eRejectHit) {
      return result;
    }
    /** Both precision and non-precision hits were tested against the same line
     * (no sync within a layer): sigma^2(non-precision test) = sigma^2(line) +
     * sigma^2(non-precision hit) sigma^2(precision test)     ~ sigma^2(line)
     * (precision hit error neglected) If the line is more precise than the
     * non-precision hit sigma(line) < sigma(non-precision hit)  <=>  2
     * sigma^2(precision test) < sigma^2(non-precision test) compatibility of
     * both hits with the line implies they belong to the same crossing: only
     * the precision hit is kept (replacing the non-precision one) without
     * branching. Otherwise the line cannot tell a wrong precision hit from the
     * right one within the non-precision hit's extent: keep both hypotheses. */
    const double precSigma{testPrec ? result.second.sigma
                                    : pat.lastHitTestRes.sigma};
    const double nonPrecSigma{testPrec ? pat.lastHitTestRes.sigma
                                       : result.second.sigma};
    if (2. * square(precSigma) < square(nonPrecSigma)) {
      ACTS_VERBOSE(__func__ << "() The line is more precise than the "
                               "non-precision hit - keep the precision hit.");
      result.first = testPrec ? LineTestDecision::eReplaceHit
                              : LineTestDecision::eRejectHit;
    } else {
      ACTS_VERBOSE(__func__ << "() The line is not more precise than the "
                               "non-precision hit - branch pattern.");
    }
    return result;
  }
  return makeResult(LineTestDecision::eBranchPattern);
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
bool GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::isPhiCompatible(
    const GeometryContext& gctx, const PatternStateAux& pat,
    const Hit_t& hit) const {
  /** We check that the test hit is compatible with the pattern phi, if
   * available. If not, we check that the test hit is in the same pattern
   * sector. We expect the test hit to be in the same hemisphere as the pattern.
   */
  const Vector3 pos{hit.globalPosition(gctx)};
  if (pat.nPhiLayers > 0) {
    if (std::abs(pat.angleToBendPlane(pos)) > m_cfg.phiTolerance) {
      ACTS_VERBOSE(__func__ << "() The pattern with phi = " << inDeg(pat.patPhi)
                            << " is not compatible with the test hit with phi "
                            << inDeg(VectorHelpers::phi(pos))
                            << " and angle to the bending plane of "
                            << inDeg(pat.angleToBendPlane(pos)));
      return false;
    }
  } else {
    if (!pat.sector.insideSector(VectorHelpers::phi(pos))) {
      ACTS_VERBOSE(__func__ << "() The test hit with phi = "
                            << inDeg(VectorHelpers::phi(pos))
                            << " is not inside the pattern sector: "
                            << static_cast<int>(pat.sector.sector()));
      return false;
    }
  }
  return true;
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
bool GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::passPatternCuts(
    const PatternStateAux& pat) const {
  /** Check that the pattern meets the minimum requirements for trigger and
   * precision layers */
  if (pat.nTriggerLayers < m_cfg.minStripEtaLayers ||
      pat.nPrecisionLayers < m_cfg.minPrecisionLayers) {
    ACTS_VERBOSE(
        __func__ << "() Pattern " << PatternStateAux::detailed(pat)
                 << "\ndoes not meet minimum layer requirements - reject.");
    return false;
  }
  /** Check requirement on the residual */
  if (pat.getMeanResidual2() > m_cfg.meanNormRes2Cut) {
    ACTS_VERBOSE(
        __func__ << "() Pattern " << PatternStateAux::detailed(pat)
                 << "\ndoes not meet the mean norm residual2 cut - reject.");
    return false;
  }
  if (countGoodGroups(pat) < m_cfg.minGroups) {
    ACTS_VERBOSE(__func__ << "() Pattern " << PatternStateAux::detailed(pat)
                          << "\ndoes not meet the minimum number of groups "
                             "with enough hits - reject.");
    return false;
  }
  return true;
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
bool GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::insertIfBest(
    const GeometryContext& gctx, PatternStateAuxVec& patterns,
    PatternStateAux&& cand) const {
  /** @brief Determine best pattern of two overlapping patterns */
  auto isBetterOverlap = [&](const PatternStateAux& a,
                             const PatternStateAux& b) {
    const int nGoodGroupDiff{countGoodGroups(a) - countGoodGroups(b)};
    if (nGoodGroupDiff != 0) {
      return nGoodGroupDiff > 0;
    }
    return isBetter(a, b);
  };
  std::vector<std::size_t> losers{};
  for (const auto& [otherIdx, other] : Acts::enumerate(patterns)) {
    if (!areOverlapping(gctx, other, cand)) {
      continue;
    }
    if (isBetterOverlap(other, cand)) {
      ACTS_VERBOSE(__func__ << "() Pattern " << PatternStateAux::detailed(other)
                            << "\nis BETTER than "
                            << PatternStateAux::detailed(cand));
      return false;
    }
    ACTS_VERBOSE(__func__ << "() Pattern " << PatternStateAux::detailed(cand)
                          << "\nis BETTER than "
                          << PatternStateAux::detailed(other));
    losers.push_back(otherIdx);
  }
  for (auto it = losers.rbegin(); it != losers.rend(); ++it) {
    if (*it + 1u != patterns.size()) {
      patterns[*it] = std::move(patterns.back());
    }
    patterns.pop_back();
  }
  patterns.push_back(std::move(cand));
  return true;
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
bool GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::areOverlapping(
    const GeometryContext& gctx, const PatternStateAux& a,
    const PatternStateAux& b) const {
  /** Check first the geometrical overlap */
  if (!a.sector.isNeighbour(b.sector)) {
    return false;
  }
  /** Check the angular difference between the seed hits */
  if (std::abs(a.patTheta - b.patTheta) > m_cfg.maxThetaSizeOverlap) {
    return false;
  }
  // Common reference radius: where the patterns actually live, e.g. the anchor
  // hit of either one.
  const double refR{
      a.lineAnchorHit->globalPosition(gctx).template head<2>().norm()};
  if (a.nPhiLayers > 1 && b.nPhiLayers > 1) {
    if (std::abs(Acts::detail::difference_periodic(
            a.phiAtRadius(refR), b.phiAtRadius(refR), 2. * std::numbers::pi)) >
        2. * m_cfg.phiTolerance) {
      return false;
    }
  } else if (a.nPhiLayers > 0) {
    if (!b.sector.insideSector(a.phiAtRadius(refR))) {
      return false;
    }
  } else if (b.nPhiLayers > 0) {
    if (!a.sector.insideSector(b.phiAtRadius(refR))) {
      return false;
    }
  }
  /** If we reach here, the patterns can overlap geometrically, so check the hit
   * content */
  std::size_t nSharedHits{0}, nSharedGroups{0}, nCommonGroups{0};
  for (std::size_t st{0u}; st < Topology_t::nGroups; ++st) {
    const auto& hitsA{a.hitsPerGroup[st]};
    const auto& hitsB{b.hitsPerGroup[st]};
    if (hitsA.empty() || hitsB.empty()) {
      continue;
    }
    ++nCommonGroups;
    const std::size_t nSharedInGroup =
        std::ranges::count_if(hitsA, [&](const OrderedHit& hitA) {
          return std::ranges::any_of(hitsB, [&hitA](const OrderedHit& hitB) {
            return hitA == hitB || hitA->overlaps(*hitB);
          });
        });
    nSharedHits += nSharedInGroup;
    /** A group is shared if the majority of the smaller pattern's hits in it
     * are shared */
    if (2u * nSharedInGroup > std::min(hitsA.size(), hitsB.size())) {
      ++nSharedGroups;
    }
  }
  /** Overlap if at least half of the hits of the smaller pattern are shared, in
   * enough shared groups */
  const std::size_t minHits{std::min(a.nBendingLayers(), b.nBendingLayers())};
  return 2u * nSharedHits >= minHits &&
         nSharedGroups >= std::min<std::size_t>(m_cfg.minGroups, nCommonGroups);
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
bool GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::isBetter(
    const PatternStateAux& a, const PatternStateAux& b) {
  const double resA{a.getMeanResidual2()};
  const double resB{b.getMeanResidual2()};
  const double resDiff{std::abs(resA - resB) / std::max(resA, resB)};
  const int nLayerDiff{a.nBendingLayers() - b.nBendingLayers()};
  const int nPrecLayDiff{a.nPrecisionLayers - b.nPrecisionLayers};

  /** For patterns that differ by 1–2 layers, don't sacrifice fit quality
   *  unless the extra layers are genuinely comparable. For a ≥3-layer
   * difference, the multiplicity advantage is strong enough to dominate. */
  if ((nLayerDiff == 0 && nPrecLayDiff == 0) ||
      (std::abs(nLayerDiff) < 2 && resDiff > 0.1 &&
       std::min(resA, resB) > 1.)) {
    return resA < resB;
  }
  if (nLayerDiff == 0) {
    return nPrecLayDiff > 0;
  }
  return nLayerDiff > 0;
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
unsigned GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::countPatterns(
    const PatternStateAuxVec& patterns, const Hit_t& seed,
    const SearchTree_t::coordinate_t& seedCoords) const {
  constexpr auto thetaIdx{Acts::toUnderlying(HitCoords::eTheta)};
  constexpr auto sectorIdx{Acts::toUnderlying(HitCoords::eSector)};
  return std::ranges::count_if(patterns, [&](const PatternStateAux& pattern) {
    if (std::abs(pattern.patTheta - seedCoords[thetaIdx]) >
            m_cfg.maxThetaSizeOverlap ||
        !pattern.sector.isNeighbour(Sector_t{
            static_cast<typename Sector_t::Index_t>(seedCoords[sectorIdx])})) {
      return false;
    }
    return pattern.isInPattern(seed);
  });
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
typename Topology_t::GroupIdx
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::countGoodGroups(
    const PatternStateAux& pat) const {
  return std::ranges::count_if(pat.hitsPerGroup, [this](const auto& hits) {
    return hits.size() >= m_cfg.minGroupLayers;
  });
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
template <OnlyPhiHitsProvider<Hit_t, Topology_t, Sector_t> OnlyPhiProvider_t>
void GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::addPhiOnlyHits(
    const GeometryContext& gctx, const BeamspotInfo& beamspotInfo,
    const OnlyPhiProvider_t& onlyPhiProvider,
    PatternStateAuxVec& patterns) const {
  PatternStateAuxVec survivingPatterns{};
  survivingPatterns.reserve(patterns.size());
  for (PatternStateAux& pat : patterns) {
    /** We look for phi-only hits in the buckets associated with the pattern */
    ACTS_VERBOSE(__func__ << "() Search for phi-only hits for pattern: "
                          << PatternStateAux::brief(pat));

    std::optional<typename Topology_t::GroupIdx> patternLineGroup{std::nullopt};

    auto projOntoPhiPlane = [&pat](const Vector3& pos) -> Vector3 {
      return pos - (pos.dot(pat.bendPlaneNorm) - pat.patPhiOffset) *
                       pat.bendPlaneNorm;
    };

    std::array<std::vector<Hit_t>, Topology_t::nGroups> phiHitsPerGroup{
        onlyPhiProvider.getPhiOnlyHits(pat, gctx)};

    for (auto [group, bucket] : Acts::enumerate(phiHitsPerGroup)) {
      for (Hit_t& newHit : bucket) {
        ACTS_VERBOSE(__func__ << "() *** Test phi-only hit "
                              << Acts::toString(*newHit.spacePoint())
                              << " in group " << static_cast<int>(group));

        /** Reject hits on a layer where an eta hit already measures phi */
        const auto& groupHits{pat.hitsPerGroup[group]};
        if (std::ranges::any_of(groupHits, [&](const OrderedHit& h) {
              return h->spacePoint()->measuresLoc0() &&
                     Topology_t::sameLayer(*h, newHit);
            })) {
          ACTS_VERBOSE(__func__ << "() Found eta hit measuring also phi on "
                                   "same layer - skip hit.");
          continue;
        }

        if (!isPhiCompatible(gctx, pat, newHit)) {
          ACTS_VERBOSE(__func__ << "() Phi-only hit not compatible");
          continue;
        }

        if (!patternLineGroup || *patternLineGroup != group) {
          if (!pat.defineLineAroundHit(gctx, newHit)) {
            ACTS_VERBOSE(__func__ << "() Invalid projection model for group "
                                  << group << " - skip hit.");
            continue;
          }
          patternLineGroup = group;
        }
        const Vector3 newHitPos{newHit.globalPosition(gctx)};
        using SeedingAux =
            Acts::Experimental::detail::CompSpacePointAuxiliaries;
        const double stripHalfLength{
            std::sqrt(newHit.spacePoint()->covariance()[Acts::toUnderlying(
                SeedingAux::ResidualIdx::bending)])};
        const Vector3 sensorDir{newHit.globalSensorDirection(gctx)};
        const Vector3 stripLow{
            projOntoPhiPlane(newHitPos - stripHalfLength * sensorDir)};
        const Vector3 stripHigh{
            projOntoPhiPlane(newHitPos + stripHalfLength * sensorDir)};
        const Vector3 stripDirOnPlane{(stripHigh - stripLow).normalized()};

        const double stripIntersect{
            Acts::detail::LineHelper::lineIntersect<3>(
                pat.linePos, pat.lineDir, stripLow, stripDirOnPlane)
                .pathLength()};
        const double stripProjLength{(stripHigh - stripLow).norm()};

        ACTS_VERBOSE(__func__ << "() Intersect distance from lower strip edge: "
                              << stripIntersect
                              << ", proj strip length: " << stripProjLength);

        /** Tolerance on the strip edges */
        constexpr double margin{10._mm};
        if (stripIntersect < -margin ||
            stripIntersect > (stripProjLength + margin)) {
          ACTS_VERBOSE(__func__ << "() The pattern falls outside the test hit "
                                   "strip in eta - skip hit.");
          continue;
        }

        /** Keep only the phi-only hit closest to the bending plane on each
         * layer */
        if (const auto it = std::ranges::find_if(pat.phiOnlyHits,
                                                 [&](const Hit_t& h) {
                                                   return Topology_t::sameLayer(
                                                       h, newHit);
                                                 });
            it != pat.phiOnlyHits.end()) {
          if (std::abs(pat.angleToBendPlane(it->globalPosition(gctx))) <
              std::abs(pat.angleToBendPlane(newHitPos))) {
            ACTS_VERBOSE(__func__ << "() A phi hit with better pull exists on "
                                     "same layer - skip hit.");
            continue;
          }
          *it = std::move(newHit);
        } else {
          pat.phiOnlyHits.push_back(std::move(newHit));
        }
        pat.updateBendPlane(gctx, beamspotInfo);
      }
    }
    if (pat.nPhiLayers < m_cfg.minPhiLayers) {
      ACTS_VERBOSE(
          __func__
          << "() Pattern " << PatternStateAux::detailed(pat) << " has only "
          << static_cast<int>(pat.nPhiLayers)
          << " phi layers, below the minimum required - reject this pattern.");
      continue;
    }
    survivingPatterns.push_back(std::move(pat));
  }
  std::swap(patterns, survivingPatterns);
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::GlobalPatternFinder(
    Config&& config, std::unique_ptr<const Logger> logger)
    : m_cfg{std::move(config)}, m_logger{std::move(logger)} {
  static_assert(std::is_move_assignable_v<PatternStateAux>);
  static_assert(std::is_move_constructible_v<PatternStateAux>);
  static_assert(std::is_copy_assignable_v<PatternStateAux>);
  static_assert(std::is_copy_constructible_v<PatternStateAux>);
  static_assert(std::is_nothrow_move_constructible_v<PatternStateAux>);
  static_assert(std::is_nothrow_move_assignable_v<PatternStateAux>);
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::OutputPattern
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::convertToPattern(
    PatternStateAux&& candidate) const {
  OutputPattern output{
      static_cast<typename Sector_t::Index_t>(candidate.sector.sector())};
  for (typename Topology_t::GroupIdx g = 0u; g < Topology_t::nGroups; ++g) {
    output.hitsPerGroup[g].reserve(candidate.hitsPerGroup[g].size());
    std::ranges::transform(candidate.hitsPerGroup[g],
                           std::back_inserter(output.hitsPerGroup[g]),
                           [](const OrderedHit& hit) { return hit.hitPtr; });
  }
  output.phiOnlyHits = std::move(candidate.phiOnlyHits);
  output.meanNormResidual2 = candidate.getMeanResidual2();
  output.patTheta = candidate.patTheta;
  output.patPhi = candidate.patPhi;
  output.sector = Sector_t{candidate.sector.sector()};
  output.nPrecisionLayers = candidate.nPrecisionLayers;
  output.nTriggerLayers = candidate.nTriggerLayers;
  output.nPhiLayers = candidate.nPhiLayers;

  return output;
}

template <GlobPatFinderHit Hit_t, SectorType Sector_t,
          PatternTopology<Hit_t> Topology_t>
std::vector<
    typename GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::OutputPattern>
GlobalPatternFinder<Hit_t, Sector_t, Topology_t>::convertToPattern(
    PatternStateAuxVec&& candidates) const {
  std::vector<OutputPattern> outPatterns{};
  outPatterns.reserve(candidates.size());

  for (PatternStateAux& pat : candidates) {
    outPatterns.push_back(convertToPattern(std::move(pat)));
  }
  return outPatterns;
}

}  // namespace Acts::Experimental::detail
