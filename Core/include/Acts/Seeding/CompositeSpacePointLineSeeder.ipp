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
    if (hitVec.size() > m_cfg.busyLayerLimit) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": The layer "
                            << layerIndex.value()
                            << " is too busy for seeding: " <<hitVec.size()
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
    SeedOptions<UncalibCont_t, CalibCont_t, Delegate_t>& options) const {
  if (!options.delegate) {
    throw std::invalid_argument(
        "CompositeSpacePointLineSeeder::nextSeed() -  Please construct a "
        "delegate");
  }

  Selector_t<UncalibCont_t> selectDelegate{};
  /// The layer hast not yet been initialized
  if (!options.m_upperLayer || !options.m_lowerLayer) {

  }

  return std::nullopt;
}

#ifdef STONJEK

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointSorter<Cont_t> Splitter_t,
          CompositeSpacePointContainer CalibCont_t,
          CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
bool CompositeSpacePointLineSeeder::prepareSeedOptions(
    SeedOptions<Cont_t, Splitter_t, CalibCont_t, Calibrator_t>& options) const {
  const auto& strawLayers{options.splitter->strawHits()};

  if (strawLayers.empty()) {
    ACTS_DEBUG(__func__ << "() " << __LINE__
                        << ": No straw hits available for seeding.");
    return false;
  }
  ACTS_DEBUG(__func__ << "():" << __LINE__ << " N straw layers "
                      << strawLayers.size() << " and  N strip layer "
                      << options.splitter->stripHits().size());

  if (std::ranges::any_of(strawLayers, [this](const Cont_t& layerHits) {
        return layerHits.size() > m_cfg.busyLayerLimit;
      })) {
    options.startWithPattern = false;
  }

  options.upperLayer = strawLayers.size() - 1;
  for (uint i_layer{0}; i_layer < strawLayers.size(); ++i_layer) {
    const Cont_t& layerHits = strawLayers[i_layer];
    ACTS_DEBUG("Layer " << i_layer << " has " << layerHits.size()
                        << " straw hits ");
  }

  while (options.lowerLayer < options.upperLayer) {
    const Cont_t& lowerLayerHits =
        options.splitter->strawHits()[options.lowerLayer];
    if (lowerLayerHits.size() > m_cfg.busyLayerLimit ||
        !firstGoodHit(lowerLayerHits, options.selector,
                      options.lowerHitIndex)) {
      ACTS_DEBUG("Skipping lower layer " << options.lowerLayer << " with "
                                         << lowerLayerHits.size() << " hits ");
      ++options.lowerLayer;
    } else {
      break;
    }
  }
  while (options.lowerLayer < options.upperLayer) {
    const Cont_t& upperLayerHits =
        options.splitter->strawHits()[options.upperLayer];
    if (upperLayerHits.size() > m_cfg.busyLayerLimit ||
        !firstGoodHit(upperLayerHits, options.selector,
                      options.upperHitIndex)) {
      --options.upperLayer;
    } else {
      break;
    }
  }
  return true;
}

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointSorter<Cont_t> Splitter_t,
          CompositeSpacePointContainer CalibCont_t,
          CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
void CompositeSpacePointLineSeeder::moveToNextCandidate(
    SeedOptions<Cont_t, Splitter_t, CalibCont_t, Calibrator_t>& options) const {
  const Cont_t& lower = options.splitter->strawHits()[options.lowerLayer];
  const Cont_t& upper = options.splitter->strawHits()[options.upperLayer];

  // Vary the left right solutions
  ++options.signComboIndex;
  if (options.signComboIndex < s_signCombo.size()) {
    return;
  }
  // All sign combos tested. Let's reset the signs combo and move on
  /// to the next hit inside the layer
  options.signComboIndex = 0;

  // Move to the next hit in the lower layer
  if (moveToNextHit(lower, options.selector, options.lowerHitIndex)) {
    return;
  }

  // Reset to the first good hit in the lower layer
  if (firstGoodHit(lower, options.selector, options.lowerHitIndex)) {
    // --> so we can update to the next hit in the upper layer
    if (moveToNextHit(upper, options.selector, options.upperHitIndex)) {
      return;
    }
  }

  // All combinations of hits & lines in both layers have been processed
  // Switch to next lower layer but skip the busy ones according to the
  // configuration
  while (options.lowerLayer < options.upperLayer) {
    const Cont_t& nextLower =
        options.splitter->strawHits()[++options.lowerLayer];
    // skip noisy layers
    if (nextLower.size() > m_cfg.busyLayerLimit) {
      continue;
    }
    // find first good hit in the new lower layer
    if (firstGoodHit(nextLower, options.selector, options.lowerHitIndex)) {
      break;
    }
  }

  // If we found a new good lower layer, reset to the first good hit in the
  // upper layer
  if (options.lowerLayer < options.upperLayer) {
    firstGoodHit(options.splitter->strawHits()[options.upperLayer],
                 options.selector, options.upperHitIndex);
    return;
  }

  // if we found a seed we would like to abort seeding after a certain lower
  // layer index, e.g. the multilayer boundary
  if (options.abortSelector.connected() &&
      options.abortSelector(options.lowerLayer) && options.nGenSeeds) {
    // Client requested to abort seeding at this layer index
    options.lowerLayer = options.upperLayer;
    return;
  }

  // Since we did not find a good seed to now, we are just going to search for
  // the first good hit in two layers as far apart as possible

  options.lowerLayer = 0;

  do {
    const Cont_t& nextLower = options.splitter->strawHits()[options.lowerLayer];
    if (nextLower.size() > m_cfg.busyLayerLimit) {
      continue;
    }
    if (firstGoodHit(nextLower, options.selector, options.lowerHitIndex)) {
      break;
    }
  } while (++options.lowerLayer < options.upperLayer);

  while (options.lowerLayer < options.upperLayer) {
    const Cont_t& nextUpper =
        options.splitter->strawHits()[--options.upperLayer];
    if (nextUpper.size() > m_cfg.busyLayerLimit) {
      continue;
    }
    if (firstGoodHit(nextUpper, options.selector, options.upperHitIndex)) {
      break;
    }
  }
}

template <typename T>
using SeedSolutionType =
    std::optional<CompositeSpacePointLineSeeder::SeedSolution<T>>;

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointSorter<Cont_t> Splitter_t,
          CompositeSpacePointContainer CalibCont_t,
          CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
SeedSolutionType<CalibCont_t> CompositeSpacePointLineSeeder::nextSeed(
    SeedOptions<Cont_t, Splitter_t, CalibCont_t, Calibrator_t>& options) const {
  std::optional<SeedSolution<Cont_t>> found{std::nullopt};
  if (!options.nGenSeeds && options.startWithPattern) {
    ++options.nGenSeeds;
    found = std::make_optional(SeedSolution<CalibCont_t>());
    found->lineParams = options.patternParams;
    found->seedHits.reserve(options.splitter->strawHits().size() +
                            options.splitter->stripHits().size());

    Vector posForCalib =
        Vector(found->lineParams[toUnderlying(ParIdx::x0)],
               found->lineParams[toUnderlying(ParIdx::y0)], 0.);
    Vector dirForCalib = makeDirectionFromPhiTheta(
        found->lineParams[toUnderlying(ParIdx::phi)],
        found->lineParams[toUnderlying(ParIdx::theta)]);
    for (const auto& strawLayerHits : options.splitter->strawHits()) {
      Cont_t tmpCalibHits = options.calibrator->calibrate(
          *options.calibContext, posForCalib, dirForCalib,
          options.patternParams[toUnderlying(ParIdx::t0)], strawLayerHits);
      found->seedHits.insert(found->seedHits.end(),
                             std::make_move_iterator(tmpCalibHits.begin()),
                             std::make_move_iterator(tmpCalibHits.end()));
    }

    found->nStrawHits = found->seedHits.size();
    for (const auto& stripLayerHits : options.splitter->stripHits()) {
      Cont_t tmpCalibHits = options.calibrator->calibrate(
          *options.calibContext, posForCalib, dirForCalib,
          options.patternParams[toUnderlying(ParIdx::t0)], stripLayerHits);
      found->seedHits.insert(found->seedHits.end(),
                             std::make_move_iterator(tmpCalibHits.begin()),
                             std::make_move_iterator(tmpCalibHits.end()));
    }

    found->y0 = options.patternParams[toUnderlying(ParIdx::y0)];
    found->theta = options.patternParams[toUnderlying(ParIdx::theta)];
    found->solutionSigns.resize(options.splitter->strawHits().size());
    options.seenSolutions.push_back(*found);

    return found;
  }
  ACTS_DEBUG("Will start looking for seeds now");
  while (options.lowerLayer < options.upperLayer) {
    found = buildSeed(options);
    moveToNextCandidate(options);
    if (found) {
      return found;
    }
  }
  return std::nullopt;
}

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
