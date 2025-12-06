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
  using ResidualIdx = detail::CompSpacePointAuxiliaries::ResidualIdx;
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
  const double distTubes = Acts::fastHypot(dY, dZ);
  assert(distTubes > 1._mm);
  constexpr auto covIdx = Acts::toUnderlying(ResidualIdx::bending);
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
  const UnCalibCont_t& strawLayer = m_splitter->strawHits().at(layIdx);
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
        const UnCalibCont_t& strawLayer = m_splitter->strawHits().at(layIdx);
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
    ostr << "    **** " << (h + 1ul) << ") " << toString(getHit(h))
         << std::endl;
  }
}

/// ##########################################################################
///                CompositeSpacePointLineSeeder::SeedOptions
/// ##########################################################################
template <
    CompositeSpacePointContainer UncalibCont_t,
    CompositeSpacePointContainer CalibCont_t,
    detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t> Delegate_t>
void CompositeSpacePointLineSeeder::SeedOptions<
    UncalibCont_t, CalibCont_t, Delegate_t>::print(std::ostream& ostr) const {
  ostr << "Seed options:\n";
  const std::size_t nStraw = delegate->strawHits().size();
  ostr << "N strawLayers: " << nStraw
       << " N strip layers: " << delegate->stripHits().size() << "\n";
  ostr << "upperLayer " << m_upperLayer.value_or(nStraw - 1ul) << " lowerLayer "
       << m_lowerLayer.value_or(0u) << " upperHitIndex " << m_upperHitIndex
       << " lower layer hit index " << m_lowerHitIndex << " sign combo index "
       << toString(encodeAmbiguity(s_signCombo[m_signComboIndex][0],
                                   s_signCombo[m_signComboIndex][1]))
       << "\n";
  ostr << " start with pattern " << startWithPattern << " nGenSeeds "
       << nGenSeeds << " nStrawCut " << nStrawCut << "\n";
}

template <CompositeSpacePointContainer UnCalibCont_t>
bool CompositeSpacePointLineSeeder::moveToNextHit(
    const UnCalibCont_t& hitVec, const Selector_t<UnCalibCont_t>& selector,
    std::size_t& hitIdx) const {
  ACTS_VERBOSE(__func__ << "() " << __LINE__
                        << ": Moving to next good hit from index " << hitIdx
                        << " in straw layer with " << hitVec.size()
                        << " hits.");
  while (hitIdx < hitVec.size()) {
    ++hitIdx;
    if (selector(*hitVec[hitIdx])) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Moved towards index "
                            << hitIdx);
      return true;
    }
  }
  return false;
}

template <CompositeSpacePointContainer UnCalibCont_t>
bool CompositeSpacePointLineSeeder::firstGoodHit(
    const UnCalibCont_t& hitVec, const Selector_t<UnCalibCont_t>& selector,
    std::size_t& hitIdx) const {
  hitIdx = 0;
  if (hitVec.empty()) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Layer is empty.");
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
  if (!layerIndex) {
    layerIndex = moveForward ? 0u : strawLayers.size() - 1u;
    const UnCalibCont_t& hitVec{strawLayers.at(layerIndex.value())};
    if (hitVec.size() <= m_cfg.busyLayerLimit &&
        firstGoodHit(hitVec, selector, hitIdx)) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Instantiated "
                            << (moveForward ? "lower" : "upper") << " layer to "
                            << layerIndex.value() << ".");
      return true;
    }
  }
  ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Move "
                        << (moveForward ? "lower" : "upper") << " layer "
                        << layerIndex.value() << " to next value.");
  /// Increment or decrement the layer index
  while ((moveForward ? (++layerIndex.value()) : (--layerIndex.value())) <
         strawLayers.size()) {
    const UnCalibCont_t& hitVec{strawLayers.at(layerIndex.value())};
    /// Check whether the layer index is still witihin the allow boundaries
    if ((moveForward && layerIndex.value() >= boundary) ||
        (!moveForward && layerIndex.value() <= boundary)) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": The "
                            << (moveForward ? "lower" : "upper") << " index "
                            << layerIndex.value()
                            << " exceeds the boundary: " << boundary << ".");
      return false;
    }
    const std::size_t nHits =
        !m_cfg.busyLimitCountGood
            ? hitVec.size()
            : std::ranges::count_if(
                  hitVec, [&](const auto& hit) { return selector(*hit); });
    if (nHits > m_cfg.busyLayerLimit) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": The layer "
                            << layerIndex.value()
                            << " is too busy for seeding: " << nHits
                            << ". Limit: " << m_cfg.busyLayerLimit << ".");
      continue;
    }

    /// Check whether a good hit can be detected inside the layer
    if (firstGoodHit(strawLayers.at(layerIndex.value()), selector, hitIdx)) {
      ACTS_VERBOSE(__func__
                   << "() " << __LINE__ << ": Loop over all hits in the  "
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
    SeedOptions<UncalibCont_t, CalibCont_t, Delegate_t>& options) const {
  if (!options.delegate) {
    throw std::invalid_argument(
        "CompositeSpacePointLineSeeder::nextSeed() -  Please construct a "
        "delegate");
  }

  const StrawLayers_t<UncalibCont_t>& strawLayers{
      options.delegate->strawHits()};

  Selector_t<UncalibCont_t> selector{};
  /// The layer hast not yet been initialized
  if (!options.m_upperLayer || !options.m_lowerLayer) {
    /// Check whether the seeding can start with the external pattern parameters
    if (options.startWithPattern &&
        std::ranges::any_of(strawLayers,
                            [this](const UncalibCont_t& layerHits) {
                              return layerHits.size() > m_cfg.busyLayerLimit;
                            })) {
      options.startWithPattern = false;
    }
    if (options.startWithPattern) {
      SegmentSeed<CalibCont_t> patternSeed{
          options.patternParams, options.delegate->newContainer(cctx)};

      const auto [pos, dir] = makeLine(options.patternParams);
      const double t0 = patternSeed.params[toUnderlying(ParIdx::t0)];

      auto append = [&](StrawLayers_t<UncalibCont_t>& hitLayers) {
        for (const auto& layer : hitLayers) {
          for (const auto& hit : layer) {
            options.delegate->append(cctx, pos, dir, t0, *hit,
                                     patternSeed.hits);
          }
        }
      };
      append(strawLayers);
      append(patternSeed->stripHits());
      options.startWithPattern = false;
      return patternSeed;
    }
    /// No valid seed can be found
    if (!nextLayer(strawLayers, selector, strawLayers.size(),
                   options.m_lowerLayer, options.m_lowerHitIndex, true) ||
        !nextLayer(strawLayers, selector, options.m_lowerHitIndex.value(),
                   options.m_upperLayer, options.m_upperHitIndex, false)) {
      ACTS_DEBUG(__func__ << "() " << __LINE__
                          << ": No valid seed can be constructed. ");
      return false;
    }
  }

  while (options.m_lowerHitIndex.value() < options.m_upperHitIndex.value()) {
    moveToNextCandidate(selector, options);
  }

  return std::nullopt;
}
template <
    CompositeSpacePointContainer UncalibCont_t,
    CompositeSpacePointContainer CalibCont_t,
    detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t> Delegate_t>
void CompositeSpacePointLineSeeder::moveToNextCandidate(
    const Selector_t<UncalibCont_t>& selector,
    SeedOptions<UncalibCont_t, CalibCont_t, Delegate_t>& options) const {
  // Vary the left right solutions
  ++options.signComboIndex;
  if (options.signComboIndex < s_signCombo.size()) {
    return;
  }
  // All sign combos tested. Let's reset the signs combo and move on
  /// to the next hit inside the layer
  options.signComboIndex = 0;
  const StrawLayers_t<UncalibCont_t>& strawLayers{
      options.delegate->strawHits()};

  const UncalibCont_t& lower = strawLayers[options.m_lowerLayer.value()];
  const UncalibCont_t& upper = strawLayers[options.m_upperLayer.value()];

  /// Next good hit in the lower layer found
  if (moveToNextHit(lower, selector, options.lowerHitIndex)) {
    return;
  }
  /// Reset the hit in the lower layer && move to the next hit
  /// in the upper layer
  if (firstGoodHit(lower, selector, options.lowerHitIndex) &&
      moveToNextHit(upper, selector, options.upperHitIndex)) {
    return;
  }
  /// All hit combinations in the two layers are tried -> move to next layer
  auto& layerToStay{options.m_moveUpLayer ? options.m_lowerLayer
                                          : options.m_upperLayer};
  auto& layerToMove{options.m_moveUpLayer ? options.m_upperLayer
                                          : options.m_lowerLayer};
  auto& hitToStay{options.m_moveUpLayer ? options.lowerHitIndex
                                        : options.upperHitIndex};
  auto& hitToMove{options.m_moveUpLayer ? options.upperHitIndex
                                        : options.lowerHitIndex};

  /// Reset the hits in the layer that remains and go to the next layer on the
  /// other side. Next time the otherside will stay and the former will move.
  if (firstGoodHit(strawLayers[layerToStay.value()], selector, hitToStay) &&
      nextLayer(strawLayers, selector, layerToStay.value(), layerToMove,
                hitToMove, options.m_moveUpLayer)) {
    options.m_moveUpLayer = !options.m_moveUpLayer;
    return;
  }
}
#ifdef STONJEK

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointSorter<Cont_t> Splitter_t,
          CompositeSpacePointContainer CalibCont_t,
          CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
SeedSolutionType<CalibCont_t> CompositeSpacePointLineSeeder::buildSeed(
    SeedOptions<Cont_t, Splitter_t, CalibCont_t, Calibrator_t>& options) const {
  auto& upperHit = options.splitter->strawHits()[options.upperLayer].at(
      options.upperHitIndex);
  auto& lowerHit = options.splitter->strawHits()[options.lowerLayer].at(
      options.lowerHitIndex);
  TangentAmbi ambi = encodeAmbiguity(s_signCombo[options.signComboIndex][0],
                                     s_signCombo[options.signComboIndex][1]);
  const CalibrationContext* ctx = options.calibContext;

  auto seedPars = constructTangentLine(*lowerHit, *upperHit, ambi);
  ACTS_DEBUG("Line pars " << seedPars << " seed opts " << options);
  if (!isValidLine(seedPars)) {
    return std::nullopt;
  }

  seedPars.lineParams =
      constructLine(seedPars.theta, seedPars.y0, options.patternParams);

  Vector posForCalib =
      Vector(seedPars.lineParams[toUnderlying(ParIdx::x0)],
             seedPars.lineParams[toUnderlying(ParIdx::y0)], 0.);
  Vector dirForCalib = makeDirectionFromPhiTheta(
      seedPars.lineParams[toUnderlying(ParIdx::phi)],
      seedPars.lineParams[toUnderlying(ParIdx::theta)]);
  if (m_cfg.recalibSeedCircles) {
    Cont_t lowerUpperToCalib{lowerHit, upperHit};
    CalibCont_t calibLowerUpper = options.calibrator->calibrate(
        *ctx, posForCalib, dirForCalib,
        options.patternParams[toUnderlying(ParIdx::t0)], lowerUpperToCalib);
    auto seedSolCalib = constructTangentLine(*(calibLowerUpper.at(0)),
                                             *(calibLowerUpper.at(1)), ambi);
    if (!isValidLine(seedSolCalib)) {
      return std::nullopt;
    }
  }
  ACTS_DEBUG("Found N seeds so far " << options.seenSolutions.size());
  // check if we have already seen this solution
  if (std::ranges::any_of(options.seenSolutions, [&seedPars](const auto& seen) {
        const double deltaY = Acts::abs(seen.y0 - seedPars.y0);
        const double limitY = Acts::fastHypot(seen.dY0, seedPars.dY0);
        const double deltaTheta = Acts::abs(seen.theta - seedPars.theta);
        const double limitTheta = Acts::fastHypot(seen.dTheta, seedPars.dTheta);
        return deltaY < limitY && deltaTheta < limitTheta;
      })) {
    return std::nullopt;
  }
  ACTS_DEBUG("start looking for compatible hits ");
  // now we search for the uncalibrated hits that are compatible with the seed
  SeedSolution<Cont_t> seedSol(seedPars);
  for (const auto [layerNr, hitsInLayer] :
       Acts::enumerate(options.splitter->strawHits())) {
    bool hadGoodHit{false};
    for (const auto& [hitNr, testMe] : Acts::enumerate(hitsInLayer)) {
      const double distance =
          Acts::abs(Acts::detail::LineHelper::signedDistance(
              testMe->localPosition(), testMe->sensorDirection(), posForCalib,
              dirForCalib));
      const double chi2 = detail::CompSpacePointAuxiliaries::chi2Term(
          posForCalib, dirForCalib, *testMe);

      ACTS_DEBUG("Hit in layer " << layerNr << " pull " << std::sqrt(chi2)
                                 << " distance " << distance << " drift radius "
                                 << testMe->driftRadius());
      if (chi2 < Acts::pow(m_cfg.hitPullCut, 2u) &&
          distance < options.strawRadius) {
        hadGoodHit = true;
        seedSol.seedHits.emplace_back(testMe);
        seedSol.nStrawHits +=
            options.selector.connected() && options.selector(*testMe);
      } else if (hadGoodHit) {
        break;
      }
    }
  }
  // check if we collected enough straw hits
  const double hitCut =
      std::max(1.0 * m_cfg.nStrawHitCut,
               m_cfg.nStrawLayHitCut * options.splitter->strawHits().size());
  ACTS_DEBUG("Found " << seedSol.nStrawHits
                      << " compatible straw hits. Hit cut is " << hitCut);
  if (seedSol.nStrawHits < hitCut) {
    return std::nullopt;
  }

  if (m_cfg.overlapCorridor) {
    seedSol.solutionSigns = detail::CompSpacePointAuxiliaries::strawSigns(
        posForCalib, dirForCalib, seedSol.seedHits);
    for (unsigned int a = options.startWithPattern;
         a < options.seenSolutions.size(); ++a) {
      const auto& acceptedSol = options.seenSolutions[a];
      unsigned int nOverlap{0};
      std::vector<int> corridor = detail::CompSpacePointAuxiliaries::strawSigns(
          posForCalib, dirForCalib, acceptedSol.seedHits);
      for (unsigned int l = 0; l < acceptedSol.seedHits.size(); ++l) {
        nOverlap += (corridor[l] == acceptedSol.solutionSigns[l]);
      }
      if (nOverlap == corridor.size() &&
          acceptedSol.seedHits.size() >= seedSol.seedHits.size()) {
        return std::nullopt;
      }
    }
  }

  // Calibrate the seed hits to be returned
  auto finalSeedSol = std::make_optional(SeedSolution<CalibCont_t>(seedPars));
  finalSeedSol->seedHits = options.calibrator->calibrate(
      *ctx, posForCalib, dirForCalib,
      options.patternParams[toUnderlying(ParIdx::t0)], seedSol.seedHits);
  finalSeedSol->nStrawHits = finalSeedSol->seedHits.size();

  // Keep track of the seed that
  options.seenSolutions.emplace_back(std::move(seedSol));
  /** If we found a long straw seed, then ensure that all
   *  subsequent seeds have at least the same amount of straw hits. */
  if (m_cfg.tightenHitCut) {
    options.nStrawCut =
        std::max(m_cfg.nStrawHitCut,
                 std::max(finalSeedSol->nStrawHits, options.nStrawCut));
  }

  ++options.nGenSeeds;

  /** Associate strip hits to the seed */

  double bestChi2Loc0{std::numeric_limits<double>::max()},
      bestChi2Loc1{std::numeric_limits<double>::max()};
  std::size_t bestIdxLoc0{0}, bestIdxLoc1{0};
  for (const auto& stripLayerHits : options.splitter->stripHits()) {
    for (const auto& [hitIdx, testMe] : Acts::enumerate(stripLayerHits)) {
      const double chi2 = detail::CompSpacePointAuxiliaries::chi2Term(
          posForCalib, dirForCalib, *testMe);
      if (testMe->measuresLoc0() && chi2 < bestChi2Loc0) {
        bestChi2Loc0 = chi2;
        bestIdxLoc0 = hitIdx;
      }
      if (testMe->measuresLoc1() && chi2 < bestChi2Loc1) {
        bestChi2Loc1 = chi2;
        bestIdxLoc1 = hitIdx;
      }
      std::vector<typename Cont_t::value_type> stripHitsToCalibrate{};
      if (bestChi2Loc0 < Acts::pow(m_cfg.hitPullCut, 2)) {
        stripHitsToCalibrate.emplace_back(stripLayerHits.at(bestIdxLoc0));
      }
      if (bestChi2Loc1 < Acts::pow(m_cfg.hitPullCut, 2) &&
          bestIdxLoc1 != bestIdxLoc0) {
        stripHitsToCalibrate.emplace_back(stripLayerHits.at(bestIdxLoc1));
      }

      auto calibratedStripHits = options.calibrator->calibrate(
          *ctx, posForCalib, dirForCalib,
          options.patternParams[toUnderlying(ParIdx::t0)],
          stripHitsToCalibrate);
      finalSeedSol->seedHits.insert(
          finalSeedSol->seedHits.end(),
          std::make_move_iterator(calibratedStripHits.begin()),
          std::make_move_iterator(calibratedStripHits.end()));
    }
  }

  return finalSeedSol;
}

#endif
}  // namespace Acts::Experimental
