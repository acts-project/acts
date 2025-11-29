// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/CompositeSpacePointLineSeeder.hpp"

#include "Acts/Surfaces/detail/LineHelper.hpp"

#include <algorithm>

namespace Acts::Experimental {

constexpr CompositeSpacePointLineSeeder::TangentAmbi
CompositeSpacePointLineSeeder::encodeAmbiguity(const int signTop,
                                               const int signBottom) {
  assert(Acts::abs(signTop) == 1 && Acts::abs(signBottom) == 1);
  using enum TangentAmbi;
  if (signTop == 1 && signBottom == 1) {
    return RR;
  } else if (signTop == 1 && signBottom == -1) {
    return RL;
  } else if (signTop == -1 && signBottom == 1) {
    return LR;
  }
  return LL;
}

constexpr CompositeSpacePointLineSeeder::TangentAmbi
CompositeSpacePointLineSeeder::encodeAmbiguity(const std::array<int, 2> ambi) {
  return encodeAmbiguity(ambi[0], ambi[1]);
}

inline std::string CompositeSpacePointLineSeeder::toString(
    const CompositeSpacePointLineSeeder::TangentAmbi ambi) {
  switch (ambi) {
    using enum TangentAmbi;
    case RR:
      return "Right - Right";
    case RL:
      return "Right - Left";
    case LR:
      return "Left - Right";
    case LL:
      return "Left - Left";
  }
  return "Undefined";
}

template <CompositeSpacePoint Spt_t>
CompositeSpacePointLineSeeder::SeedParameters
CompositeSpacePointLineSeeder::constructTangentLine(const Spt_t& topHit,
                                                    const Spt_t& bottomHit,
                                                    const TangentAmbi ambi) {
  using ResidualIdx = detail::CompSpacePointAuxiliaries::ResidualIdx;
  SeedParameters result{};
  const auto& [signTop, signBot] = s_signCombo[toUnderlying(ambi)];

  const Vector& bottomPos{bottomHit.localPosition()};
  const Vector& topPos{topHit.localPosition()};
  const Vector& eY{bottomHit.toNextSensor()};
  const Vector& eZ{bottomHit.planeNormal()};
  const Vector D = topPos - bottomPos;

  assert(Acts::abs(eY.dot(eZ)) < std::numeric_limits<double>::epsilon());
  assert(Acts::abs(bottomHit.sensorDirection().dot(eY)) <
         std::numeric_limits<double>::epsilon());
  assert(Acts::abs(bottomHit.sensorDirection().dot(eZ)) <
         std::numeric_limits<double>::epsilon());

  assert(topHit.isStraw() && bottomHit.isStraw());
  const double dY = D.dot(eY);
  const double dZ = D.dot(eZ);

  const double thetaTubes = std::atan2(dY, dZ);
  const double distTubes = Acts::fastHypot(dY, dZ);

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
  if (denomSquare < std::numeric_limits<double>::epsilon()) {
    return result;
  }
  result.dTheta = combDriftUncert / std::sqrt(denomSquare) / distTubes;
  result.dY0 =
      std::hypot(bottomPos.dot(eY) * sinTheta + bottomPos.dot(eZ) * cosTheta,
                 1.) *
      result.dTheta;
  if(result.theta < 0){
    result.theta += std::numbers::pi;
  }
  return result;
}

template <CompositeSpacePointContainer UnCalibCont_t>
bool CompositeSpacePointLineSeeder::moveToNextHit(
    const UnCalibCont_t& hitVec,
    const Selector_t<SpacePoint_t<UnCalibCont_t>>& selector,
    std::size_t& hitIdx) const {
  ACTS_DEBUG(__func__ << "() " << __LINE__
                         << ": Moving to next good hit from index " << hitIdx
                         << " in hit vector of size " << hitVec.size()
                        << " with selector " << selector.connected());
  if(hitIdx+1 < hitVec.size()){
    ACTS_DEBUG(" Next hit is good "<< !selector(*hitVec[hitIdx+1]));
  }
  while (++hitIdx < hitVec.size() && selector.connected() &&
         !selector(*hitVec[hitIdx])) {
  }
  return hitIdx < hitVec.size() && selector.connected() &&
         selector(*hitVec[hitIdx]);
}

template <CompositeSpacePointContainer UnCalibCont_t>
bool CompositeSpacePointLineSeeder::firstGoodHit(
    const UnCalibCont_t& hitVec,
    const Selector_t<SpacePoint_t<UnCalibCont_t>>& selector,
    std::size_t& hitIdx) const {
  hitIdx = 0;
  if (hitVec.empty()) return false;
  return !selector.connected() || selector(*hitVec[hitIdx]) || moveToNextHit(hitVec, selector, hitIdx);
}

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointSorter<Cont_t> Splitter_t,
          CompositeSpacePointContainer CalibCont_t,
          CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
bool CompositeSpacePointLineSeeder::prepareSeedOptions(
    SeedOptions<Cont_t, Splitter_t, CalibCont_t, Calibrator_t>& options) const {
  if (options.splitter->strawHits().empty()) {
    ACTS_DEBUG(__func__ << "() " << __LINE__
                        << ": No straw hits available for seeding.");
    return false;
  }
  ACTS_DEBUG(__func__ << "():" << __LINE__ << " N straw layers "
                      << options.splitter->strawHits().size()
                      << " and  N strip layer "
                      << options.splitter->stripHits().size());

  if (std::ranges::any_of(options.splitter->strawHits(),
                          [this](const Cont_t& layerHits) {
                            return layerHits.size() > m_cfg.busyLayerLimit;
                          })) {
    options.startWithPattern = false;
  }

  options.upperLayer = options.splitter->strawHits().size() - 1;
  for (uint i_layer{0}; i_layer < options.splitter->strawHits().size();
       ++i_layer) {
    const Cont_t& layerHits = options.splitter->strawHits()[i_layer];
    ACTS_DEBUG("Layer " << i_layer << " has " << layerHits.size()
                        << " straw hits ");
  }

  while (options.lowerLayer < options.upperLayer) {
    const Cont_t& lowerLayerHits =
        options.splitter->strawHits()[options.lowerLayer];
    if (lowerLayerHits.size() > m_cfg.busyLayerLimit ||
        !firstGoodHit(lowerLayerHits, options.selector,
                      options.lowerHitIndex)) {
      ACTS_DEBUG("Skipping lower layer " << options.lowerLayer
                                           << " with "
                                           << lowerLayerHits.size()
                                           << " hits ");
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
  if (++options.m_signComboIndex < s_signCombo.size()) {
    return;
  }
  // All sign combos tested. Let's reset the signs combo and move on to the next
  // layer
  options.m_signComboIndex = 0;

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
};

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
    found->line = Line_t(options.patternParams);
    found->seedHits.reserve(options.splitter->strawHits().size() +
                            options.splitter->stripHits().size());

    for (const auto& strawLayerHits : options.splitter->strawHits()) {
      Cont_t tmpCalibHits = options.calibrator->calibrate(
          *options.calibContext, found->line.position(),
          found->line.direction(), options.t0Estimate, strawLayerHits);
      found->seedHits.insert(found->seedHits.end(),
                             std::make_move_iterator(tmpCalibHits.begin()),
                             std::make_move_iterator(tmpCalibHits.end()));
    }

    found->nStrawHits = found->seedHits.size();
    for (const auto& stripLayerHits : options.splitter->stripHits()) {
      Cont_t tmpCalibHits = options.calibrator->calibrate(
          *options.calibContext, found->line.position(),
          found->line.direction(), options.t0Estimate, stripLayerHits);
      found->seedHits.insert(found->seedHits.end(),
                             std::make_move_iterator(tmpCalibHits.begin()),
                             std::make_move_iterator(tmpCalibHits.end()));
    }

    found->y0 = options.patternParams[toUnderlying(Line_t::ParIndex::y0)];
    found->theta = options.patternParams[toUnderlying(Line_t::ParIndex::theta)];
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
  TangentAmbi ambi = encodeAmbiguity(s_signCombo[options.m_signComboIndex]);
  const CalibrationContext* ctx = options.calibContext;

  auto seedPars = constructTangentLine(*lowerHit, *upperHit, ambi);
  ACTS_DEBUG("Line pars " <<  seedPars << " seed opts "<< options);
  if (!isValidLine(seedPars)) {
    return std::nullopt;
  }

  seedPars.line =
      constructLine(seedPars.theta, seedPars.y0, options.patternParams);

  if (m_cfg.recalibSeedCircles) {
    Cont_t lowerUpperToCalib{lowerHit, upperHit};
    CalibCont_t calibLowerUpper = options.calibrator->calibrate(
        *ctx, seedPars.line.position(), seedPars.line.direction(),
        options.t0Estimate, lowerUpperToCalib);
    auto seedSolCalib = constructTangentLine(*(calibLowerUpper.at(0)),
                                             *(calibLowerUpper.at(1)), ambi);
    if (!isValidLine(seedSolCalib)) {
      return std::nullopt;
    }
  }
  ACTS_DEBUG("Found N seeds so far " << options.seenSolutions.size());
  // check if we have already seen this solution
  if (std::ranges::find_if(
          options.seenSolutions, [&seedPars](const auto& seen) {
            const double deltaY = Acts::abs(seen.y0 - seedPars.y0);
            const double limitY = std::hypot(seen.dY0, seedPars.dY0);
            const double deltaTheta = Acts::abs(seen.theta - seedPars.theta);
            const double limitTheta = std::hypot(seen.dTheta, seedPars.dTheta);
            return deltaY < limitY && deltaTheta < limitTheta;
          }) != options.seenSolutions.end()) {
    return std::nullopt;
  }
  ACTS_DEBUG("start looking for compatible hits ");
  // now we search for the uncalibrated hits that are compatible with the seed
  SeedSolution<Cont_t> seedSol(seedPars);
  for (const auto [layerNr, hitsInLayer] :
       Acts::enumerate(options.splitter->strawHits())) {
    bool hadGoodHit{false};
    for (const auto& testMe : hitsInLayer) {
      const double distance =
          Acts::abs(Acts::detail::LineHelper::signedDistance(
              testMe->localPosition(), testMe->sensorDirection(),
              seedSol.line.position(), seedSol.line.direction()));
      const double pull =
          Acts::abs(distance - testMe->driftRadius()) /
          std::sqrt(testMe->covariance()[Acts::toUnderlying(
              detail::CompSpacePointAuxiliaries::ResidualIdx::bending)]);

      if (pull < m_cfg.hitPullCut && distance < testMe->driftRadius()) {
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
  const unsigned hitCut =
      std::max(1. * m_cfg.nStrawHitCut,
               m_cfg.nStrawLayHitCut * options.splitter->strawHits().size());
  ACTS_DEBUG("Found " << seedSol.nStrawHits << " compatible straw hits. Hit cut is " << hitCut);
  if (seedSol.nStrawHits < hitCut) {
    return std::nullopt;
  }

  if (m_cfg.overlapCorridor) {
    seedSol.solutionSigns = detail::CompSpacePointAuxiliaries::strawSigns(
        seedSol.line, seedSol.seedHits);
    for (unsigned int a = options.startWithPattern;
         a < options.seenSolutions.size(); ++a) {
      const auto& acceptedSol = options.seenSolutions[a];
      unsigned int nOverlap{0};
      std::vector<int> corridor = detail::CompSpacePointAuxiliaries::strawSigns(
          seedSol.line, acceptedSol.seedHits);
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
      *ctx, seedSol.line.position(), seedSol.line.direction(),
      options.t0Estimate, seedSol.seedHits);
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
  return finalSeedSol;

  /** Associate strip hits */
}

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointSorter<Cont_t> Splitter_t,
          CompositeSpacePointContainer CalibCont_t,
          CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
std::ostream& CompositeSpacePointLineSeeder::SeedOptions<
    Cont_t, Splitter_t, CalibCont_t, Calibrator_t>::print(std::ostream& ostr)
    const {
  ostr << "Seed options:\n";
  ostr << "N strawLayers: " << splitter->strawHits().size()
       << " N strip layers: " << splitter->stripHits().size() << "\n";
  ostr << "upperLayer " << upperLayer << " lowerLayer " << lowerLayer
       << " upperHitIndex " << upperHitIndex << " lower layer hit index "
       << lowerHitIndex << " sign combo index " << m_signComboIndex << "\n";
  ostr << " start with pattern " << startWithPattern << " nGenSeeds "
       << nGenSeeds << " nStrawCut " << nStrawCut << "\n";
  return ostr;
}

}  // namespace Acts::Experimental
