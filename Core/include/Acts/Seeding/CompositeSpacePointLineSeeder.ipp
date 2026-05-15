// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/CompositeSpacePointLineSeeder.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <cassert>

namespace Acts::Experimental {

template <CompositeSpacePoint Sp_t>
CompositeSpacePointLineSeeder::TwoCircleTangentPars
CompositeSpacePointLineSeeder::constructTangentLine(const Sp_t& topHit,
                                                    const Sp_t& bottomHit,
                                                    const TangentAmbi ambi) {
  using namespace Acts::UnitLiterals;
  using namespace Acts::detail;
  TwoCircleTangentPars result{};
  result.ambi = ambi;
  const auto& [signTop, signBot] = s_signCombo[toUnderlying(ambi)];

  const Vector& bottomPos{bottomHit.localPosition()};
  const Vector& topPos{topHit.localPosition()};
  const Vector& eY{bottomHit.toNextSensor()};
  const Vector& eZ{bottomHit.planeNormal()};
  const Vector D = topPos - bottomPos;

  assert(Acts::abs(eY.dot(eZ)) < s_epsilon);
  assert(Acts::abs(bottomHit.sensorDirection().dot(eY)) < s_epsilon);
  assert(Acts::abs(bottomHit.sensorDirection().dot(eZ)) < s_epsilon);
  assert(topHit.isStraw() && bottomHit.isStraw());

  const double dY = D.dot(eY);
  const double dZ = D.dot(eZ);

  const double thetaTubes = std::atan2(dY, dZ);
  const double distTubes = fastHypot(dY, dZ);
  assert(distTubes > 1._mm);
  constexpr auto covIdx = toUnderlying(CovIdx::bending);
  const double combDriftUncert{topHit.covariance()[covIdx] +
                               bottomHit.covariance()[covIdx]};
  const double R =
      -signBot * bottomHit.driftRadius() + signTop * topHit.driftRadius();
  result.theta = thetaTubes - std::asin(std::clamp(R / distTubes, -1., 1.));

  const double cosTheta = std::cos(result.theta);
  const double sinTheta = std::sin(result.theta);

  result.y0 = bottomPos.dot(eY) * cosTheta - bottomPos.dot(eZ) * sinTheta -
              signBot * bottomHit.driftRadius();
  assert(Acts::abs(topPos.dot(eY) * cosTheta - topPos.dot(eZ) * sinTheta -
                   signTop * topHit.driftRadius() - result.y0) <
         std::numeric_limits<float>::epsilon());
  result.y0 /= cosTheta;
  const double denomSquare = 1. - Acts::pow(R / distTubes, 2);
  if (denomSquare < s_epsilon) {
    return result;
  }
  result.dTheta = combDriftUncert / std::sqrt(denomSquare) / distTubes;
  result.dY0 =
      Acts::fastHypot(
          bottomPos.dot(eY) * sinTheta + bottomPos.dot(eZ) * cosTheta, 1.) *
      result.dTheta;

  return result;
}

template <CompositeSpacePoint Sp_t>
CompositeSpacePointLineSeeder::Vector
CompositeSpacePointLineSeeder::makeDirection(const Sp_t& refHit,
                                             const double tanAngle) {
  const Vector& eY{refHit.toNextSensor()};
  const Vector& eZ{refHit.planeNormal()};
  const double cosTheta = std::cos(tanAngle);
  const double sinTheta = std::sin(tanAngle);
  return copySign<Vector, double>(sinTheta * eY + cosTheta * eZ, sinTheta);
}

template <CompositeSpacePointContainer Cont_t>
std::size_t CompositeSpacePointLineSeeder::countHits(
    const Cont_t& container, const Selector_t<Cont_t>& selector) const {
  if (m_cfg.busyLimitCountGood) {
    return std::ranges::count_if(
        container, [&](const auto& hit) { return selector(*hit); });
  }
  return container.size();
}

/// ##########################################################################
///                CompositeSpacePointLineSeeder::SeedSolution
/// ##########################################################################
template <CompositeSpacePointContainer UnCalibCont_t,
          detail::CompositeSpacePointSorter<UnCalibCont_t> Splitter_t>
CompositeSpacePointLineSeeder::SeedSolution<UnCalibCont_t,
                                            Splitter_t>::SpacePoint_t
CompositeSpacePointLineSeeder::SeedSolution<UnCalibCont_t, Splitter_t>::getHit(
    const std::size_t idx) const {
  const auto& [layIdx, hitIdx] = m_seedHits.at(idx);
  const UnCalibCont_t& strawLayer = m_splitter.strawHits().at(layIdx);
  return *strawLayer.at(hitIdx);
}
template <CompositeSpacePointContainer UnCalibCont_t,
          detail::CompositeSpacePointSorter<UnCalibCont_t> Splitter_t>
std::vector<int> CompositeSpacePointLineSeeder::SeedSolution<
    UnCalibCont_t, Splitter_t>::leftRightAmbiguity(const Vector& seedPos,
                                                   const Vector& seedDir)
    const {
  std::vector<int> result{};
  result.reserve(size());
  std::ranges::transform(
      m_seedHits, std::back_inserter(result),
      [&](const std::pair<std::size_t, std::size_t>& indices) {
        const auto& [layIdx, hitIdx] = indices;
        const UnCalibCont_t& strawLayer = m_splitter.strawHits().at(layIdx);
        return detail::CompSpacePointAuxiliaries::strawSign(
            seedPos, seedDir, *strawLayer.at(hitIdx));
      });
  return result;
}
template <CompositeSpacePointContainer UnCalibCont_t,
          detail::CompositeSpacePointSorter<UnCalibCont_t> Splitter_t>
void CompositeSpacePointLineSeeder::SeedSolution<
    UnCalibCont_t, Splitter_t>::append(const std::size_t layIdx,
                                       const std::size_t hitIdx) {
  assert(!rangeContainsValue(m_seedHits, std::make_pair(layIdx, hitIdx)));
  m_seedHits.emplace_back(layIdx, hitIdx);
}
template <CompositeSpacePointContainer UnCalibCont_t,
          detail::CompositeSpacePointSorter<UnCalibCont_t> Splitter_t>
void CompositeSpacePointLineSeeder::SeedSolution<
    UnCalibCont_t, Splitter_t>::print(std::ostream& ostr) const {
  TwoCircleTangentPars::print(ostr);
  const std::size_t N = size();
  ostr << ", associated hits: " << N << std::endl;
  for (std::size_t h = 0; h < N; ++h) {
    ostr << "    **** " << (h + 1ul) << ") " << Acts::toString(getHit(h))
         << std::endl;
  }
}

/// ##########################################################################
///                CompositeSpacePointLineSeeder::SeedingState
/// ##########################################################################
template <
    CompositeSpacePointContainer UncalibCont_t,
    CompositeSpacePointContainer CalibCont_t,
    detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t> Delegate_t>
void CompositeSpacePointLineSeeder::SeedingState<
    UncalibCont_t, CalibCont_t, Delegate_t>::print(std::ostream& ostr) const {
  const std::size_t nStraw = this->strawHits().size();

  ostr << "Seed state:\n";
  ostr << "N strawLayers: " << nStraw
       << " N strip layers: " << this->stripHits().size() << "\n";
  ostr << "upperLayer " << m_upperLayer.value_or(nStraw - 1ul) << " lowerLayer "
       << m_lowerLayer.value_or(0u) << " upperHitIndex " << m_upperHitIndex
       << " lower layer hit index " << m_lowerHitIndex << " sign combo index "
       << toString(encodeAmbiguity(s_signCombo[m_signComboIndex][0],
                                   s_signCombo[m_signComboIndex][1]))
       << "\n";
  ostr << "Number of seeds " << nGenSeeds() << " nStrawCut " << m_nStrawCut
       << "\n";
  if (nGenSeeds() > 0ul) {
    for (const auto& seen : m_seenSolutions) {
      ostr << "################################################" << std::endl;
      ostr << seen << std::endl;
      ostr << "################################################" << std::endl;
    }
  }
}

template <CompositeSpacePointContainer UnCalibCont_t>
bool CompositeSpacePointLineSeeder::moveToNextHit(
    const UnCalibCont_t& hitVec, const Selector_t<UnCalibCont_t>& selector,
    std::size_t& hitIdx) const {
  ACTS_VERBOSE(__func__ << "() " << __LINE__
                        << " - Moving to next good hit from index " << hitIdx
                        << " in straw layer with " << hitVec.size()
                        << " hits.");
  while (++hitIdx < hitVec.size()) {
    if (selector(*hitVec[hitIdx])) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Moved towards index "
                            << hitIdx);
      return true;
    }
  }
  ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - No good hit found.");
  return false;
}

template <CompositeSpacePointContainer UnCalibCont_t>
bool CompositeSpacePointLineSeeder::firstGoodHit(
    const UnCalibCont_t& hitVec, const Selector_t<UnCalibCont_t>& selector,
    std::size_t& hitIdx) const {
  hitIdx = 0;
  if (hitVec.empty()) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Layer is empty.");
    return false;
  }
  return selector(*hitVec[hitIdx]) || moveToNextHit(hitVec, selector, hitIdx);
}
template <CompositeSpacePointContainer UnCalibCont_t>
bool CompositeSpacePointLineSeeder::nextLayer(
    const StrawLayers_t<UnCalibCont_t>& strawLayers,
    const Selector_t<UnCalibCont_t>& selector, const std::size_t boundary,
    std::optional<std::size_t>& layerIndex, std::size_t& hitIdx,
    bool moveForward) const {
  if (strawLayers.empty()) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - No straw layers.");
    return false;
  }
  /// The layer index is not yet instantiated.
  if (!layerIndex.has_value()) {
    layerIndex = moveForward ? 0u : strawLayers.size() - 1u;
    const UnCalibCont_t& hitVec{strawLayers.at(layerIndex.value())};
    if (hitVec.size() <= m_cfg.busyLayerLimit &&
        firstGoodHit(hitVec, selector, hitIdx)) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Instantiated "
                            << (moveForward ? "lower" : "upper") << " layer to "
                            << layerIndex.value() << ".");
      return true;
    }
  }
  ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Move "
                        << (moveForward ? "lower" : "upper") << " layer "
                        << layerIndex.value() << " to next value.");
  /// Increment or decrement the layer index
  while ((moveForward ? (++layerIndex.value()) : (--layerIndex.value())) <
         strawLayers.size()) {
    const UnCalibCont_t& hitVec{strawLayers.at(layerIndex.value())};
    /// Check whether the layer index is still witihin the allow boundaries
    if ((moveForward && layerIndex.value() >= boundary) ||
        (!moveForward && layerIndex.value() <= boundary)) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - The "
                            << (moveForward ? "lower" : "upper") << " index "
                            << layerIndex.value()
                            << " exceeds the boundary: " << boundary << ".");
      return false;
    }

    if (const std::size_t nHits = countHits(hitVec, selector);
        nHits > m_cfg.busyLayerLimit) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - The layer "
                            << layerIndex.value()
                            << " is too busy for seeding: " << nHits
                            << ". Limit: " << m_cfg.busyLayerLimit << ".");
      continue;
    }

    /// Check whether a good hit can be detected inside the layer
    if (firstGoodHit(strawLayers.at(layerIndex.value()), selector, hitIdx)) {
      ACTS_VERBOSE(__func__
                   << "() " << __LINE__ << " - Loop over all hits in the  "
                   << layerIndex.value() << (moveForward ? "lower" : "upper")
                   << " layer.");
      return true;
    }
  }
  return false;
}

template <
    CompositeSpacePointContainer UncalibCont_t,
    CompositeSpacePointContainer CalibCont_t,
    detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t> Delegate_t>
std::optional<CompositeSpacePointLineSeeder::SegmentSeed<CalibCont_t>>
CompositeSpacePointLineSeeder::nextSeed(
    const CalibrationContext& cctx,
    SeedingState<UncalibCont_t, CalibCont_t, Delegate_t>& state) const {
  const StrawLayers_t<UncalibCont_t>& strawLayers{state.strawHits()};

  Selector_t<UncalibCont_t> selector{};
  selector.template connect<&Delegate_t::goodCandidate>(&state);
  /// The layer hast not yet been initialized
  if (!state.m_upperLayer || !state.m_lowerLayer) {
    state.m_nStrawCut = m_cfg.nStrawHitCut;
    /// Check whether the seeding can start with the external pattern
    /// parameters
    if (!state.m_patternSeedProduced &&
        (!m_cfg.startWithPattern ||
         std::ranges::any_of(strawLayers, [&](const UncalibCont_t& layerHits) {
           return countHits(layerHits, selector) > m_cfg.busyLayerLimit;
         }))) {
      state.m_patternSeedProduced = true;
    }
    if (!state.m_patternSeedProduced) {
      SegmentSeed<CalibCont_t> patternSeed{state.initialParameters(),
                                           state.newContainer(cctx)};

      const auto [pos, dir] = makeLine(state.initialParameters());
      const double t0 = patternSeed.parameters[toUnderlying(ParIdx::t0)];

      auto append = [&](const StrawLayers_t<UncalibCont_t>& hitLayers) {
        for (const auto& layer : hitLayers) {
          for (const auto& hit : layer) {
            state.append(cctx, pos, dir, t0, *hit, patternSeed.hits);
          }
        }
      };
      append(strawLayers);
      append(state.stripHits());
      state.m_patternSeedProduced = true;
      return patternSeed;
    }
    ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Instantiate layers. ");
    /// No valid seed can be found
    if (!nextLayer(strawLayers, selector, strawLayers.size(),
                   state.m_lowerLayer, state.m_lowerHitIndex, true) ||
        !nextLayer(strawLayers, selector, state.m_lowerLayer.value(),
                   state.m_upperLayer, state.m_upperHitIndex, false)) {
      ACTS_DEBUG(__func__ << "() " << __LINE__
                          << " - No valid seed can be constructed. ");
      return std::nullopt;
    }
  }

  while (state.m_lowerLayer.value() < state.m_upperLayer.value()) {
    auto seed = buildSeed(cctx, selector, state);
    moveToNextCandidate(selector, state);
    if (seed) {
      return seed;
    }
  }

  return std::nullopt;
}
template <
    CompositeSpacePointContainer UncalibCont_t,
    CompositeSpacePointContainer CalibCont_t,
    detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t> Delegate_t>
void CompositeSpacePointLineSeeder::moveToNextCandidate(
    const Selector_t<UncalibCont_t>& selector,
    SeedingState<UncalibCont_t, CalibCont_t, Delegate_t>& state) const {
  // Vary the left right solutions
  ++state.m_signComboIndex;
  if (state.m_signComboIndex < s_signCombo.size()) {
    return;
  }
  // All sign combos tested. Let's reset the signs combo and move on
  /// to the next hit inside the layer
  state.m_signComboIndex = 0;

  const StrawLayers_t<UncalibCont_t>& strawLayers{state.strawHits()};

  const UncalibCont_t& lower = strawLayers[state.m_lowerLayer.value()];
  const UncalibCont_t& upper = strawLayers[state.m_upperLayer.value()];

  /// Next good hit in the lower layer found
  ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Move to next lower hit.");
  if (moveToNextHit(lower, selector, state.m_lowerHitIndex)) {
    return;
  }

  /// Reset the hit in the lower layer && move to the next hit
  /// in the upper layer
  ACTS_VERBOSE(__func__ << "() " << __LINE__
                        << " - All lower hits were tried increment upper hit.");
  if (firstGoodHit(lower, selector, state.m_lowerHitIndex) &&
      moveToNextHit(upper, selector, state.m_upperHitIndex)) {
    return;
  }
  /// All hit combinations in the two layers are tried -> move to next layer
  auto& layerToStay{state.m_moveUpLayer ? state.m_lowerLayer
                                        : state.m_upperLayer};
  auto& layerToMove{state.m_moveUpLayer ? state.m_upperLayer
                                        : state.m_lowerLayer};
  auto& hitToStay{state.m_moveUpLayer ? state.m_lowerHitIndex
                                      : state.m_upperHitIndex};
  auto& hitToMove{state.m_moveUpLayer ? state.m_upperHitIndex
                                      : state.m_lowerHitIndex};
  ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Move towards the next "
                        << (state.m_moveUpLayer ? "upper" : "lower")
                        << " layer. Current state: " << layerToMove.value());
  /// Reset the hits in the layer that remains and go to the next layer on the
  /// other side. Next time the otherside will stay and the former will move.
  if (firstGoodHit(strawLayers[layerToStay.value()], selector, hitToStay) &&
      nextLayer(strawLayers, selector, layerToStay.value(), layerToMove,
                hitToMove, !state.m_moveUpLayer)) {
    state.m_moveUpLayer = !state.m_moveUpLayer;
    if (state.stopSeeding(state.m_lowerLayer.value(),
                          state.m_upperLayer.value())) {
      state.m_lowerLayer.value() = state.m_upperLayer.value() + 1ul;
    }
    return;
  }
}

template <
    CompositeSpacePointContainer UncalibCont_t,
    CompositeSpacePointContainer CalibCont_t,
    detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t> Delegate_t>
bool CompositeSpacePointLineSeeder::passSeedCuts(
    const Line_t& tangentSeed,
    SeedSolution<UncalibCont_t, Delegate_t>& newSolution,
    SeedingState<UncalibCont_t, CalibCont_t, Delegate_t>& state) const {
  // check if we collected enough straw hits
  const auto& [seedPos, seedDir] = tangentSeed;
  const double hitCut =
      std::max(1.0 * state.m_nStrawCut,
               m_cfg.nStrawLayHitCut * state.strawHits().size());
  ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Found "
                        << newSolution.nStrawHits
                        << " compatible straw hits. Hit cut is " << hitCut);
  if (newSolution.nStrawHits < hitCut) {
    return false;
  }
  if (!m_cfg.overlapCorridor) {
    return true;
  }
  newSolution.solutionSigns = newSolution.leftRightAmbiguity(seedPos, seedDir);
  for (std::size_t a = 0ul; a < state.nGenSeeds(); ++a) {
    const auto& acceptedSol = state.m_seenSolutions[a];
    std::size_t nOverlap{0};
    const std::vector<int> corridor =
        acceptedSol.leftRightAmbiguity(seedPos, seedDir);
    for (std::size_t l = 0; l < acceptedSol.size(); ++l) {
      nOverlap += (corridor[l] == acceptedSol.solutionSigns[l]);
    }
    if (nOverlap == corridor.size() &&
        acceptedSol.size() >= newSolution.size()) {
      return false;
    }
  }
  if (m_cfg.tightenHitCut) {
    state.m_nStrawCut = std::max(state.m_nStrawCut, newSolution.nStrawHits);
  }
  return true;
}

template <
    CompositeSpacePointContainer UncalibCont_t,
    CompositeSpacePointContainer CalibCont_t,
    detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t> Delegate_t>
std::optional<CompositeSpacePointLineSeeder::SegmentSeed<CalibCont_t>>
CompositeSpacePointLineSeeder::buildSeed(
    const CalibrationContext& cctx, const Selector_t<UncalibCont_t>& selector,
    SeedingState<UncalibCont_t, CalibCont_t, Delegate_t>& state) const {
  const StrawLayers_t<UncalibCont_t>& strawLayers{state.strawHits()};

  ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Try to draw new seed from \n"
                      << state << ".");
  const auto& upperHit =
      *strawLayers.at(state.m_upperLayer.value()).at(state.m_upperHitIndex);
  const auto& lowerHit =
      *strawLayers.at(state.m_lowerLayer.value()).at(state.m_lowerHitIndex);
  const auto ambi{static_cast<TangentAmbi>(state.m_signComboIndex)};
  ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - " << toString(ambi)
                        << "\n   Top seed hit: " << Acts::toString(upperHit)
                        << "\n   Bottom seed hit:" << Acts::toString(lowerHit)
                        << ".");

  const TwoCircleTangentPars seedPars =
      constructTangentLine(lowerHit, upperHit, ambi);
  ACTS_VERBOSE(__func__ << "() " << __LINE__
                        << " - Tangential parameters: " << seedPars);
  if (!isValidLine(seedPars)) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Reject seed.");
    return std::nullopt;
  }
  // check if we have already seen this solution
  if (std::ranges::any_of(state.m_seenSolutions, [&seedPars](const auto& seen) {
        const double deltaY = Acts::abs(seen.y0 - seedPars.y0);
        const double limitY = Acts::fastHypot(seen.dY0, seedPars.dY0);
        const double deltaTheta = Acts::abs(seen.theta - seedPars.theta);
        const double limitTheta = Acts::fastHypot(seen.dTheta, seedPars.dTheta);
        return deltaY < limitY && deltaTheta < limitTheta;
      })) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Reject seed.");
    return std::nullopt;
  }
  /// Continue to construct a new solution
  const double t0 = state.initialParameters()[toUnderlying(ParIdx::t0)];
  SeedSolution<UncalibCont_t, Delegate_t> newSolution{seedPars, state};

  ACTS_DEBUG(__func__ << "() " << __LINE__
                      << " - Start looking for compatible hits");

  const Line_t tangentSeed{seedPars.y0 * Vector3::UnitY(),
                           makeDirection(lowerHit, seedPars.theta)};
  const auto& [seedPos, seedDir] = tangentSeed;
  const double maxPullSq{Acts::square(m_cfg.hitPullCut)};
  constexpr auto covIdx = Acts::toUnderlying(CovIdx::bending);

  for (const auto& [layerNr, hitsInLayer] :
       Acts::enumerate(state.strawHits())) {
    bool hadGoodHit{false};
    for (const auto& [hitNr, testMe] : Acts::enumerate(hitsInLayer)) {
      using namespace Acts::detail::LineHelper;
      const double distance = Acts::abs(
          signedDistance(testMe->localPosition(), testMe->sensorDirection(),
                         seedPos, seedDir));
      // the hits are ordered in the layer so we assume that once we found good
      // hits we are moving away from the seed line so we can abort the hit
      // association
      const double rMax = state.strawRadius(*testMe);
      assert(rMax > Acts::s_epsilon);
      assert(testMe->covariance()[covIdx] > Acts::s_epsilon);
      if (distance < rMax) {
        const double pullSq =
            m_cfg.useSimpleStrawPull
                ? Acts::square(distance - Acts::abs(testMe->driftRadius())) /
                      testMe->covariance()[covIdx]
                : state.candidateChi2(cctx, seedPos, seedDir, t0, *testMe);

        if (pullSq < maxPullSq) {
          ACTS_VERBOSE(__func__
                       << "() " << __LINE__ << " - layer,hit = (" << layerNr
                       << "," << hitNr
                       << ") -> spacePoint: " << Acts::toString(*testMe)
                       << " is close enough: " << distance
                       << ", pull: " << std::sqrt(pullSq));
          hadGoodHit = true;
          newSolution.append(layerNr, hitNr);
          newSolution.nStrawHits += selector(*testMe);
          continue;
        }
      }
      if (hadGoodHit) {
        break;
      }
    }
  }
  if (!passSeedCuts(tangentSeed, newSolution, state)) {
    return std::nullopt;
  }
  return consructSegmentSeed(cctx, tangentSeed, state, std::move(newSolution));
}

template <
    CompositeSpacePointContainer UncalibCont_t,
    CompositeSpacePointContainer CalibCont_t,
    detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t> Delegate_t>
CompositeSpacePointLineSeeder::SegmentSeed<CalibCont_t>
CompositeSpacePointLineSeeder::consructSegmentSeed(
    const CalibrationContext& cctx, const Line_t& tangentSeed,
    SeedingState<UncalibCont_t, CalibCont_t, Delegate_t>& state,
    SeedSolution<UncalibCont_t, Delegate_t>&& newSolution) const {
  ACTS_DEBUG(__func__ << "() " << __LINE__ << " Construct new seed from \n"
                      << newSolution);
  SegmentSeed<CalibCont_t> finalSeed{
      combineWithPattern(tangentSeed, state.initialParameters()),
      state.newContainer(cctx)};
  const auto [seedPos, seedDir] = makeLine(finalSeed.parameters);
  /// Append the collected straws to the seed
  const double t0 = finalSeed.parameters[toUnderlying(ParIdx::t0)];
  for (std::size_t s = 0; s < newSolution.size(); ++s) {
    state.append(cctx, seedPos, seedDir, t0, newSolution.getHit(s),
                 finalSeed.hits);
  }
  /// The solution object is no longer needed for this seed.
  state.m_seenSolutions.push_back(std::move(newSolution));

  ACTS_DEBUG(__func__ << "() " << __LINE__
                      << " - Associate the strip hits to the seed");
  for (const auto& stripLayerHits : state.stripHits()) {
    double bestChi2Loc0{m_cfg.hitPullCut};
    double bestChi2Loc1{m_cfg.hitPullCut};
    std::size_t bestIdxLoc0{stripLayerHits.size()};
    std::size_t bestIdxLoc1{stripLayerHits.size()};
    /// Find the hit with the lowest pull
    for (const auto& [hitIdx, testMe] : Acts::enumerate(stripLayerHits)) {
      const double chi2 =
          state.candidateChi2(cctx, seedPos, seedDir, t0, *testMe);
      if (testMe->measuresLoc0() && chi2 < bestChi2Loc0) {
        bestChi2Loc0 = chi2;
        bestIdxLoc0 = hitIdx;
      }
      if (testMe->measuresLoc1() && chi2 < bestChi2Loc1) {
        bestChi2Loc1 = chi2;
        bestIdxLoc1 = hitIdx;
      }
    }
    if (bestIdxLoc0 < stripLayerHits.size()) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Append loc0 strip hit "
                            << Acts::toString(*stripLayerHits.at(bestIdxLoc0))
                            << ".");
      state.append(cctx, seedPos, seedDir, t0, *stripLayerHits.at(bestIdxLoc0),
                   finalSeed.hits);
    }
    if (bestIdxLoc1 != bestIdxLoc0 && bestIdxLoc1 < stripLayerHits.size()) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Append loc1 strip hit "
                            << Acts::toString(*stripLayerHits.at(bestIdxLoc1))
                            << ".");
      state.append(cctx, seedPos, seedDir, t0, *stripLayerHits.at(bestIdxLoc1),
                   finalSeed.hits);
    }
  }
  return finalSeed;
}
}  // namespace Acts::Experimental
