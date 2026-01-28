// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once
#include "Acts/Seeding/StrawChamberLineSeeder.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <format>

namespace Acts {

template <CompositeSpacePointContainer UnCalibCont_t,
          CompositeSpacePointSorter<UnCalibCont_t> Splitter_t,
          CompositeSpacePointContainer CalibSp_t,
          CompositeSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
StrawChamberLineSeeder<UnCalibCont_t, Splitter_t, CalibSp_t, Calibrator_t>::
    StrawChamberLineSeeder(const SeedParam_t& grainSeedPars,
                           const UnCalibCont_t& seedHits, const Config& cfg,
                           std::unique_ptr<const Logger> logObj)
    : m_coarsePars{grainSeedPars},
      m_hitLayers{seedHits},
      m_cfg{cfg},
      m_logger{std::move(logObj)} {
  /** The StrawChamberLine needs to have at least 2 straw layers */
  if (m_hitLayers.strawHits().size() < 2) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__
                          << " - Not enough straw layers have been parsed");
    return;
  }
  if (std::ranges::find_if(m_hitLayers.strawHits(),
                           [this](const UnCalibCont_t& vec) {
                             return vec.size() > m_cfg.busyLayerLimit;
                           }) != m_hitLayers.strawHits().end()) {
    ACTS_VERBOSE(
        __func__ << "() " << __LINE__
                 << " - Detected at least one busy layer with more than "
                 << m_cfg.busyLayerLimit
                 << ". Will not attempt the external seed as first seed");
    m_cfg.startWithPattern = false;
  }
  // Set the start for the upper layer
  m_upperLayer = m_hitLayers.strawHits().size() - 1;

  /** Check whether the first layer is too busy */
  while (m_lowerLayer < m_upperLayer &&
         m_hitLayers.strawHits()[m_lowerLayer].size() > m_cfg.busyLayerLimit) {
    ACTS_VERBOSE(
        __func__ << "() " << __LINE__ << " - Lower layer " << m_lowerLayer
                 << " has too many hits. Don't use it as seeding layer");
    ++m_lowerLayer;
  }
  /** Check whether the lower layer is too busy */
  while (m_lowerLayer < m_upperLayer &&
         m_hitLayers.strawHits()[m_upperLayer].size() > m_cfg.busyLayerLimit) {
    ACTS_VERBOSE(
        __func__ << "() " << __LINE__ << " - Upper layer " << m_upperLayer
                 << " has too many hits. Don't use it as seeding layer");
    --m_upperLayer;
  }
}
template <CompositeSpacePointContainer UnCalibCont_t,
          CompositeSpacePointSorter<UnCalibCont_t> Splitter_t,
          CompositeSpacePointContainer CalibSp_t,
          CompositeSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
std::ostream& StrawChamberLineSeeder<
    UnCalibCont_t, Splitter_t, CalibSp_t,
    Calibrator_t>::SeedSolution::print(std::ostream& ostr) const {
  ostr << "theta: " << theta << " pm " << dTheta << ", ";
  ostr << "y0: " << Y0 << " pm " << dY0 << ", ";
  ostr << "signs: ";
  for (std::size_t i = 0; i < solutionSigns.size(); ++i) {
    ostr << solutionSigns[i];
    if (i + 1 != solutionSigns.size()) {
      ostr << ", ";
    }
  }
  return ostr;
}
template <CompositeSpacePointContainer UnCalibCont_t,
          CompositeSpacePointSorter<UnCalibCont_t> Splitter_t,
          CompositeSpacePointContainer CalibSp_t,
          CompositeSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
std::optional<typename StrawChamberLineSeeder<
    UnCalibCont_t, Splitter_t, CalibSp_t, Calibrator_t>::DriftCircleSeed>
StrawChamberLineSeeder<UnCalibCont_t, Splitter_t, CalibSp_t, Calibrator_t>::
    generateSeed(const CalibrationContext& ctx) {
  if (!numGenerated() && m_cfg.startWithPattern) {
    DriftCircleSeed patternSeed{};
    patternSeed.parameters = m_coarsePars;
    const Vector3 seedPos{m_coarsePars[eBoundLoc0], m_coarsePars[eBoundLoc1],
                          0.};
    const Vector3 seedDir = makeDirectionFromPhiTheta(
        m_coarsePars[eBoundPhi], m_coarsePars[eBoundTheta]);

    for (const auto& hitLayers :
         {m_hitLayers.strawHits(), m_hitLayers.stripHits()}) {
      for (const auto& hits : hitLayers) {
        for (const UncalibSp_t& hit : hits) {
          m_cfg.calibrator->calibrate(ctx, seedPos, seedDir,
                                      m_coarsePars[eBoundTime], *hit,
                                      patternSeed.measurements);
        }
      }
    }
    ++m_nGenSeeds;
    return patternSeed;
  }
  std::optional<DriftCircleSeed> found{std::nullopt};

  while (m_lowerLayer < m_upperLayer) {
    const UnCalibCont_t& lower = m_hitLayers.strawHits().at(m_lowerLayer);
    const UnCalibCont_t& upper = m_hitLayers.strawHits().at(m_upperLayer);
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << " - Layers with hits: "
                          << m_hitLayers.strawHits().size()
                          << " -- next bottom hit: " << m_lowerLayer
                          << ", hit: " << m_lowerHitIndex << " ("
                          << lower.size() << "), topHit " << m_upperLayer
                          << ", " << m_upperHitIndex << " (" << upper.size()
                          << ") - ambiguity ("
                          << s_signCombos[m_signComboIndex][0] << ";"
                          << s_signCombos[m_signComboIndex][1] << ")");

    found = buildSeed(ctx, upper.at(m_upperHitIndex), lower.at(m_lowerHitIndex),
                      s_signCombos.at(m_signComboIndex));
    /// Increment for the next candidate
    moveToNextCandidate();
    /// If a candidate is built return it. Otherwise continue the process
    if (found) {
      return found;
    }
  }
  return std::nullopt;
}
template <CompositeSpacePointContainer UnCalibCont_t,
          CompositeSpacePointSorter<UnCalibCont_t> Splitter_t,
          CompositeSpacePointContainer CalibSp_t,
          CompositeSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
void StrawChamberLineSeeder<UnCalibCont_t, Splitter_t, CalibSp_t,
                            Calibrator_t>::moveToNextCandidate() {
  const UnCalibCont_t& lower = m_hitLayers.strawHits()[m_lowerLayer];
  const UnCalibCont_t& upper = m_hitLayers.strawHits()[m_upperLayer];
  /// Vary the left-right solutions
  if (++m_signComboIndex < s_signCombos.size()) {
    return;
  }
  m_signComboIndex = 0;

  /// Move to the next hit in the lower layer
  if (++m_lowerHitIndex < lower.size()) {
    return;
  }
  m_lowerHitIndex = 0;
  /// Move to the next hit in the upper layer
  if (++m_upperHitIndex < upper.size()) {
    return;
  }
  m_upperHitIndex = 0;
  /// All combinations of hits & lines in both layers are processed
  /// Switch to the next lowerLayer. But skip the busy ones according to the
  /// configuration
  while (m_lowerLayer < m_upperLayer &&
         m_hitLayers.strawHits()[++m_lowerLayer].size() >
             m_cfg.busyLayerLimit) {
  }

  if (m_lowerLayer < m_upperLayer) {
    return;
  }
  m_lowerLayer = 0;
  while (m_lowerLayer < m_upperLayer &&
         m_hitLayers.strawHits()[--m_upperLayer].size() >
             m_cfg.busyLayerLimit) {
  }
}

template <CompositeSpacePointContainer UnCalibCont_t,
          CompositeSpacePointSorter<UnCalibCont_t> Splitter_t,
          CompositeSpacePointContainer CalibSp_t,
          CompositeSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
std::optional<typename StrawChamberLineSeeder<
    UnCalibCont_t, Splitter_t, CalibSp_t, Calibrator_t>::DriftCircleSeed>
StrawChamberLineSeeder<UnCalibCont_t, Splitter_t, CalibSp_t,
                       Calibrator_t>::buildSeed(const CalibrationContext& ctx,
                                                const UncalibSp_t& topHit,
                                                const UncalibSp_t& bottomHit,
                                                const SignCombo_t& signs) {
  // Fetch the signs
  const auto& [signTop, signBot] = signs;
  /// Calculate the relative radius
  double R =
      signBot * bottomHit->driftRadius() - signTop * topHit->driftRadius();
  const Vector3& bottomPos{bottomHit->localPosition()};
  const Vector3& topPos{topHit->localPosition()};
  /// Calculate the distance between the two and their relative angle
  const Vector3 D = topPos - bottomPos;
  const double thetaTubes = std::atan2(D.y(), D.z());
  const double distTubes = fastHypot(D.y(), D.z());
  ACTS_VERBOSE(__func__ << "() " << __LINE__
                        << " - Try to build new 2 circle seed from bottom Hit: "
                        << toString(bottomPos)
                        << ", r: " << bottomHit->driftRadius() << ", top hit: "
                        << toString(topPos) << ", r: " << topHit->driftRadius()
                        << " --> tube distance: " << toString(D)
                        << ", mag: " << distTubes << ", theta: " << thetaTubes);

  /// Calculate the seed theta.
  double theta{thetaTubes - std::asin(std::clamp(R / distTubes, -1., 1.))};

  Vector3 seedDir =
      makeDirectionFromPhiTheta(90. * UnitConstants::degree, theta);
  double y0 = bottomPos.y() * seedDir.z() - bottomPos.z() * seedDir.y() +
              signBot * bottomHit->driftRadius();

  DriftCircleSeed candidateSeed{};

  double combDriftUncert{
      std::sqrt(bottomHit->covariance()[eY] + topHit->covariance()[eY])};

  candidateSeed.parameters[eBoundTime] = m_coarsePars[eBoundTime];
  candidateSeed.parameters[eBoundLoc1] = y0 / seedDir.z();

  /// Check that the initial estimate of the seed is in range
  if (!detail::inRange(m_cfg.thetaRange, theta) ||
      !detail::inRange(m_cfg.interceptRange,
                       candidateSeed.parameters[eBoundLoc1])) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__
                          << " Seed parameters are out of range");
    return std::nullopt;
  }

  Vector3 seedPos = candidateSeed.parameters[eBoundLoc1] * Vector3::UnitY();

  assert(std::abs(topPos.y() * seedDir.z() - topPos.z() * seedDir.y() +
                  signTop * topHit->driftRadius() - y0) <
         std::numeric_limits<float>::epsilon());
  ACTS_VERBOSE(__func__ << "() " << __LINE__
                        << ": Candidate seed theta: " << theta
                        << ", tanTheta: " << (seedDir.y() / seedDir.z())
                        << ", y0: " << candidateSeed.parameters[eBoundLoc1]);

  SeedSolution solCandidate{};
  solCandidate.Y0 = candidateSeed.parameters[eBoundLoc1];
  solCandidate.theta = theta;
  /** Reserve enough memory */
  solCandidate.seedHits.reserve(2 * m_hitLayers.strawHits().size());
  /// d/dx asin(x) = 1 / sqrt(1- x*x)
  const double denomSquare = 1. - std::pow(R / distTubes, 2);
  if (denomSquare < std::numeric_limits<double>::epsilon()) {
    ACTS_VERBOSE("Invalid seed, rejecting");
    return std::nullopt;
  }
  solCandidate.dTheta = combDriftUncert / std::sqrt(denomSquare) / distTubes;
  solCandidate.dY0 =
      fastHypot(-bottomPos.y() * seedDir.y() + bottomPos.z() * seedDir.z(),
                1.) *
      solCandidate.dTheta;
  ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Test new " << solCandidate
                        << ". " << m_seenSolutions.size());
  if (std::ranges::find_if(
          m_seenSolutions, [&solCandidate, this](const SeedSolution& seen) {
            const double deltaY = std::abs(seen.Y0 - solCandidate.Y0);
            const double limitY = fastHypot(seen.dY0, solCandidate.dY0);
            const double dTheta = std::abs(seen.theta - solCandidate.theta);
            const double limitTh = fastHypot(seen.dTheta, solCandidate.dTheta);
            ACTS_VERBOSE(__func__
                         << "() " << __LINE__ << ": " << seen << " delta Y: "
                         << deltaY << " " << (deltaY < limitY ? '<' : '>')
                         << " " << limitY << " delta theta: " << dTheta << " "
                         << (dTheta < limitTh ? '<' : '>') << " " << limitTh);
            return deltaY < limitY && dTheta < limitTh;
            ;
          }) != m_seenSolutions.end()) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Reject due to similarity");
    return std::nullopt;
  }
  /** Collect all hits close to the seed line */
  for (const auto& [layerNr, hitsInLayer] :
       enumerate(m_hitLayers.strawHits())) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": " << hitsInLayer.size()
                          << " hits in layer " << (layerNr + 1));
    bool hadGoodHit{false};
    for (const UncalibSp_t& testMe : hitsInLayer) {
      const double pull =
          0.;  // detail::calcStrawResidual(seedPos, seedDir, *testMe);
      ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Test hit: "
                            << toString(testMe->localPosition()) << "radius: "
                            << testMe->driftRadius() << ", pull: " << pull);
      if (pull < m_cfg.hitPullCut) {
        hadGoodHit = true;
        solCandidate.seedHits.emplace_back(testMe);
        ++candidateSeed.nStrawHits;
      }  /// what ever comes after will be further away from the segment
      else if (hadGoodHit) {
        break;
      }
    }
  }
  /** Reject seeds with too little straw hit association */
  const unsigned hitCut =
      std::max(1. * m_cfg.nStrawHitCut,
               m_cfg.nStrawLayHitCut * m_hitLayers.strawHits().size());
  if (1. * candidateSeed.nStrawHits < hitCut) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Too few hits associated "
                          << candidateSeed.nStrawHits << ", expect at least "
                          << hitCut << " hits.");
    return std::nullopt;
  }
  /* Reject seeds that are in the same corridor as previously built seed and
   * don't add any other straw hit*/
  if (m_cfg.overlapCorridor) {
    // solCandidate.solutionSigns =
    //     detail::calcStrawSigns(seedPos, seedDir, solCandidate.seedHits);
    ACTS_VERBOSE(__func__ << "() " << __LINE__
                          << ": Circle solutions for seed - " << solCandidate);
    /** Last check whether another seed with the same left-right combination
     * hasn't already been found */
    for (const SeedSolution& accepted : m_seenSolutions) {
      unsigned int nOverlap{0};
      std::vector<int> corridor = {};
      // detail::calcStrawSigns(seedPos, seedDir, accepted.seedHits);
      for (unsigned int l = 0; l < accepted.seedHits.size(); ++l) {
        nOverlap += (corridor[l] == accepted.solutionSigns[l]);
      }
      /// The seed basically generates a new line that's in he same left-right
      /// corridor compared to a previously found solution. There's no need to
      /// return that seed again
      if (nOverlap == corridor.size() &&
          accepted.seedHits.size() >= solCandidate.seedHits.size()) {
        ACTS_VERBOSE(
            __func__
            << "() " << __LINE__
            << ": Same set of hits collected within the same corridor.");
        return std::nullopt;
      }
    }
  }

  /** If we found a long straw hit seed, then ensure that all
   *  subsequent seeds have at least the same amount of straw hits hits. */
  if (m_cfg.tightenHitCut) {
    m_cfg.nStrawHitCut = std::max(m_cfg.nStrawHitCut, candidateSeed.nStrawHits);
  }
  /** Increment the counter of the total number of generated seeds */
  ++m_nGenSeeds;

  /** Combine the theta estimate from the y-z plane with the phi estimate
   *  from the external seed */
  combineWithExternalSeed(candidateSeed.parameters[eBoundLoc1], theta,
                          candidateSeed.parameters);

  /** Update the direction & position vectors to take the external parameters
   *  along the non-bending plane into account */
  seedDir = makeDirectionFromPhiTheta(candidateSeed.parameters[eBoundPhi],
                                      candidateSeed.parameters[eBoundTheta]);
  seedPos[eX] = candidateSeed.parameters[eBoundLoc0];
  using unwrapped_ref = std::unwrap_reference_t<UncalibSp_t>;
  for (unwrapped_ref onSeed : solCandidate.seedHits) {
    m_cfg.calibrator->calibrate(ctx, seedPos, seedDir,
                                candidateSeed.parameters[eBoundTime], *onSeed,
                                candidateSeed.measurements);
  }


    seedDir = makeDirectionFromPhiTheta(candidateSeed.parameters[eBoundPhi],
                                        candidateSeed.parameters[eBoundTheta]);
    seedPos[eY] = candidateSeed.parameters[eBoundLoc1];
  }

  /// Finally complement the seed with the set of best strip layer measurements
  for (const auto& [layerNr, hitsInLayer] :
       enumerate(m_hitLayers.stripHits())) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": " << hitsInLayer.size()
                          << " hits in layer " << (layerNr + 1));

    double bestRes{m_cfg.hitPullCut};
    std::size_t bestIdx = hitsInLayer.size();
    for (const auto& [hitIdx, testMe] : enumerate(hitsInLayer)) {
      const double res =
          0.;  // detail::calcStripResidual(seedPos, seedDir, *testMe);
      if (res < bestRes) {
        bestRes = res;
        bestIdx = hitIdx;
      }
    }
    /// Add the best hit to the seed
    if (bestIdx < hitsInLayer.size()) {
      m_cfg.calibrator->calibrate(
          ctx, seedPos, seedDir, candidateSeed.parameters[eBoundTime],
          *hitsInLayer[bestIdx], candidateSeed.measurements);
    }
  }
  return candidateSeed;
}

template <CompositeSpacePointContainer UnCalibCont_t,
          CompositeSpacePointSorter<UnCalibCont_t> Splitter_t,
          CompositeSpacePointContainer CalibSp_t,
          CompositeSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
void StrawChamberLineSeeder<
    UnCalibCont_t, Splitter_t, CalibSp_t,
    Calibrator_t>::combineWithExternalSeed(const double y0, const double theta,
                                           SeedParam_t& outPars) const {
  const double tanBeta = std::tan(theta);
  const double tanAlpha = std::cos(m_coarsePars[eBoundPhi]) * tanBeta;
  outPars.template block<2, 1>(eBoundPhi, 0) = makePhiThetaFromDirection(
      makeDirectionFromAxisTangents(tanAlpha, tanBeta));
  outPars[eBoundLoc0] = m_coarsePars[eBoundLoc0];
  outPars[eBoundTime] = m_coarsePars[eBoundTime];
  outPars[eBoundLoc1] = y0;
}
}  // namespace Acts
