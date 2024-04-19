// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <unordered_map>

namespace Acts {

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
const TrackContainer<track_container_t, traj_t, holder_t>
Acts::AthenaAmbiguityResolution::prepareOutputTrack(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<std::size_t>& goodTracks) const {
  auto trackStateContainer = tracks.trackStateContainerHolder();
  auto trackContainer = std::make_shared<VectorTrackContainer>();
  trackContainer->reserve(goodTracks.size());
  // Temporary empty track state container: we don't change the original one,
  // but we need one for filtering
  auto tempTrackStateContainer = std::make_shared<VectorMultiTrajectory>();
  TrackContainer solvedTracks{trackContainer, tempTrackStateContainer};
  solvedTracks.ensureDynamicColumns(tracks);

  for (auto&& iTrack : goodTracks) {
    auto destProxy = solvedTracks.getTrack(solvedTracks.addTrack());
    auto srcProxy = tracks.getTrack(iTrack);
    destProxy.copyFrom(srcProxy, false);
    destProxy.tipIndex() = srcProxy.tipIndex();
  }

  const TrackContainer<track_container_t, traj_t, holder_t> outputTracks{
      std::make_shared<VectorTrackContainer>(std::move(*trackContainer)),
      trackStateContainer};
  return outputTracks;
}

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t, typename source_link_hash_t,
          typename source_link_equality_t>
std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
AthenaAmbiguityResolution::computeInitialState(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    source_link_hash_t&& sourceLinkHash,
    source_link_equality_t&& sourceLinkEquality) const {
  auto measurementIndexMap =
      std::unordered_map<SourceLink, std::size_t, source_link_hash_t,
                         source_link_equality_t>(0, sourceLinkHash,
                                                 sourceLinkEquality);

  std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
      measurementsPerTrack;

  ACTS_VERBOSE("Starting to compute initial state");

  for (const auto& track : tracks) {
    std::vector<std::tuple<std::size_t, std::size_t, bool>> measurements_tuples;

    for (auto ts : track.trackStatesReversed()) {
      ACTS_DEBUG("Track state type check: ");
      if (ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        Acts::SourceLink sourceLink = ts.getUncalibratedSourceLink();
        ACTS_DEBUG("Track state type is MeasurementFlag");
        const auto& geoID = ts.referenceSurface().geometryId();

        // assign a new measurement index if the source link was not seen yet
        auto emplace = measurementIndexMap.try_emplace(
            sourceLink, measurementIndexMap.size());

        bool isoutliner =
            ts.typeFlags().test(Acts::TrackStateFlag::OutlierFlag);

        measurements_tuples.push_back(
            std::make_tuple(emplace.first->second, geoID.volume(), isoutliner));
      }
    }

    measurementsPerTrack.push_back(std::move(measurements_tuples));
  }

  return measurementsPerTrack;
}

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
std::vector<double> Acts::AthenaAmbiguityResolution::simpleScore(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<std::map<std::size_t, TrackFeatures>>& trackFeaturesMaps,
    Optional_cuts<track_container_t, traj_t, holder_t> optionalCuts) const {
  std::vector<double> trackScore;

  int iTrack = 0;

  ACTS_VERBOSE("Number of detectors: " << m_cfg.detectorMap.size());

  ACTS_INFO("Starting to score tracks");

  // Loop over all the tracks in the container
  for (const auto& track : tracks) {
    auto trackFeaturesMap = std::map<std::size_t, TrackFeatures>();

    // Loop over all the track states in a track for counting
    // hits/hole/outliers per detector.
    for (const auto& ts : track.trackStatesReversed()) {
      auto iTypeFlags = ts.typeFlags();

      if (iTypeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
        auto iVolume = ts.referenceSurface().geometryId().volume();
        auto volume_it = m_cfg.volumeMap.find(iVolume);
        if (volume_it == m_cfg.volumeMap.end()) {
          ACTS_ERROR("Volume " << iVolume << "not found in the volume map");
          continue;
        }
        auto detectorId = volume_it->second;
        if (iTypeFlags.test(Acts::TrackStateFlag::SharedHitFlag)) {
          trackFeaturesMap[detectorId].nSharedHits++;
        }
        trackFeaturesMap[detectorId].nHits++;
      } else if (iTypeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
        auto iVolume = ts.referenceSurface().geometryId().volume();
        auto volume_it = m_cfg.volumeMap.find(iVolume);
        if (volume_it == m_cfg.volumeMap.end()) {
          ACTS_ERROR("Volume " << iVolume << "not found in the volume map");
          continue;
        }
        auto detectorId = volume_it->second;
        trackFeaturesMap[detectorId].nHoles++;
      } else if (iTypeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
        auto iVolume = ts.referenceSurface().geometryId().volume();
        auto volume_it = m_cfg.volumeMap.find(iVolume);
        if (volume_it == m_cfg.volumeMap.end()) {
          ACTS_ERROR("Volume " << iVolume << "not found in the volume map");
          continue;
        }
        auto detectorId = volume_it->second;
        trackFeaturesMap[detectorId].nOutliers++;
      }
    }
    trackFeaturesMaps.push_back(trackFeaturesMap);

    double score = 1;
    // cuts on pT
    if (Acts::VectorHelpers::perp(track.momentum()) > m_cfg.pTMax ||
        Acts::VectorHelpers::perp(track.momentum()) < m_cfg.pTMin) {
      score = 0;
      iTrack++;
      trackScore.push_back(score);
      ACTS_DEBUG("Track: " << iTrack
                           << " has score = 0, due to pT cuts --- pT = "
                           << Acts::VectorHelpers::perp(track.momentum()));
      continue;
    }

    // cuts on phi
    if (Acts::VectorHelpers::phi(track.momentum()) > m_cfg.phiMax ||
        Acts::VectorHelpers::phi(track.momentum()) < m_cfg.phiMin) {
      score = 0;
      iTrack++;
      trackScore.push_back(score);
      ACTS_DEBUG("Track: " << iTrack
                           << " has score = 0, due to phi cuts --- phi =  "
                           << Acts::VectorHelpers::phi(track.momentum()));
      continue;
    }

    // cuts on eta
    if (Acts::VectorHelpers::eta(track.momentum()) > m_cfg.etaMax ||
        Acts::VectorHelpers::eta(track.momentum()) < m_cfg.etaMin) {
      score = 0;
      iTrack++;
      trackScore.push_back(score);
      ACTS_DEBUG("Track: " << iTrack
                           << " has score = 0, due to eta cuts --- eta =  "
                           << Acts::VectorHelpers::eta(track.momentum()));
      continue;
    }

    // cuts on optional cuts
    for (const auto& ambicut : optionalCuts.cuts) {
      if (ambicut(track)) {
        score = 0;
        ACTS_DEBUG("Track: " << iTrack
                             << " has score = 0, due to optional cuts.");
        break;
      }
    }

    if (score == 0) {
      iTrack++;
      trackScore.push_back(score);
      ACTS_DEBUG("Track: " << iTrack << " score : " << score);
      continue;
    }
    // Reject tracks which didn't pass the detector cuts.
    for (std::size_t detectorId = 0; detectorId < m_cfg.detectorMap.size();
         detectorId++) {
      auto detector_it = m_cfg.detectorMap.find(detectorId);
      auto detector = detector_it->second;

      ACTS_DEBUG("---> Found summary information");
      ACTS_DEBUG("---> Detector ID: " << detectorId);
      ACTS_DEBUG("---> Number of hits: " << trackFeaturesMap[detectorId].nHits);
      ACTS_DEBUG(
          "---> Number of holes: " << trackFeaturesMap[detectorId].nHoles);
      ACTS_DEBUG("---> Number of outliers: "
                 << trackFeaturesMap[detectorId].nOutliers);

      if ((trackFeaturesMap[detectorId].nHits < detector.minHits) ||
          (trackFeaturesMap[detectorId].nHits > detector.maxHits) ||
          (trackFeaturesMap[detectorId].nHoles > detector.maxHoles) ||
          (trackFeaturesMap[detectorId].nOutliers > detector.maxOutliers)) {
        score = 0;
        ACTS_DEBUG("Track: " << iTrack
                             << " has score = 0, due to detector cuts");
        break;
      }
    }

    if (score == 0) {
      iTrack++;
      trackScore.push_back(score);
      ACTS_DEBUG("Track: " << iTrack << " score : " << score);
      continue;
    }

    // real scoring starts here

    if (!m_cfg.useAmbiguityFunction) {
      // if the ambiguity function is used, the score is processed with a
      // different algorithm than the simple score.
      score = 100;
      // Adding the score for each detector.
      // detector score is determined by the number of hits/hole/outliers *
      // hit/hole/outlier scoreWeights in a detector.
      for (std::size_t detectorId = 0; detectorId < m_cfg.detectorMap.size();
           detectorId++) {
        auto detector_it = m_cfg.detectorMap.find(detectorId);
        auto detector = detector_it->second;
        score += trackFeaturesMap[detectorId].nHits * detector.hitsScoreWeight;
        score +=
            trackFeaturesMap[detectorId].nHoles * detector.holesScoreWeight;
        score += trackFeaturesMap[detectorId].nOutliers *
                 detector.outliersScoreWeight;
        score += trackFeaturesMap[detectorId].nSharedHits *
                 detector.otherScoreWeight;
      }

      // Adding scores based on optional weights
      for (const auto& ambiweights : optionalCuts.weights) {
        ambiweights(track, score);
      }

      // Adding the score based on the chi2/ndf
      if (track.chi2() > 0 && track.nDoF() > 0) {
        double p = 1. / log10(10. + 10. * track.chi2() / track.nDoF());
        if (p > 0) {
          score += p;
        } else
          score -= 50;
      }
    }

    iTrack++;

    // Add the score to the vector
    trackScore.push_back(score);
    ACTS_VERBOSE("Track: " << iTrack << " score: " << score);

  }  // end of loop over tracks

  if (!m_cfg.useAmbiguityFunction) {
    ACTS_VERBOSE("Not using ambiguity function");
    return trackScore;
  }

  ACTS_VERBOSE("Using ambiguity function");

  std::vector<double> trackScoreAmbig;
  iTrack = 0;
  // Loop over all the tracks in the container
  for (const auto& track : tracks) {
    double score = trackScore[iTrack];
    // remove tracks which did not pass the previous cuts.
    if (score == 0) {
      trackScoreAmbig.push_back(0.0f);
      iTrack++;
      continue;
    }
    auto trackFeaturesMap =
        trackFeaturesMaps[iTrack];  // get the trackFeatures map for the track
    double pT = Acts::VectorHelpers::perp(track.momentum());
    // start with larger score for tracks with higher pT.
    double prob = log10(pT / UnitConstants::MeV) - 1.;
    // pT in GeV, hence 100 MeV is minimum and gets score = 1
    ACTS_DEBUG("Modifier for pT = " << pT << " GeV is : " << prob
                                    << "  New score now: " << prob);

    for (std::size_t detectorId = 0; detectorId < m_cfg.detectorMap.size();
         detectorId++) {
      auto detector_it = m_cfg.detectorMap.find(detectorId);
      auto detector = detector_it->second;

      // choosing a scaling factor based on the number of hits in a track per
      // detector.
      std::size_t iHits = trackFeaturesMap[detectorId].nHits;
      if (detector.factorHits.size() < iHits) {
        ACTS_WARNING("Detector " << detectorId
                                 << " has not enough factorhits in the "
                                    "detector.factorHits vector");
        continue;
      }
      if (iHits > detector.maxHits) {
        prob = prob * (detector.maxHits - iHits + 1);  // hits are good !
        iHits = detector.maxHits;
      }
      prob = prob * detector.factorHits[iHits];
      ACTS_DEBUG("Modifier for " << iHits
                                 << " hits: " << detector.factorHits[iHits]
                                 << "  New score now: " << prob);

      // choosing a scaling factor based on the number of holes in a track per
      // detector.
      std::size_t iHoles = trackFeaturesMap[detectorId].nHoles;
      if (detector.factorHoles.size() < iHoles) {
        ACTS_WARNING("Detector " << detectorId
                                 << " has not enough factorholes in the "
                                    "detector.factorHoles vector");
        continue;
      }
      if (iHoles > detector.maxHoles) {
        prob /= (iHoles - detector.maxHoles + 1);  // holes are bad !
        iHoles = detector.maxHoles;
      }
      prob = prob * detector.factorHoles[iHoles];
      ACTS_DEBUG("Modifier for " << iHoles
                                 << " holes: " << detector.factorHoles[iHoles]
                                 << "  New score now: " << prob);
    }

    for (const auto& ambiscore : optionalCuts.ambiscores) {
      ambiscore(track, prob);
    }

    if (track.chi2() > 0 && track.nDoF() > 0) {
      double chi2 = track.chi2();
      int indf = track.nDoF();
      double fac = 1. / log10(10. + 10. * chi2 / indf);
      prob = prob * fac;
      ACTS_DEBUG("Modifier for chi2 = " << chi2 << " and NDF = " << indf
                                        << " is : " << fac
                                        << "  New score now: " << prob)
    }
    trackScoreAmbig.push_back(prob);
  }
  return trackScoreAmbig;
}

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
std::vector<int> Acts::AthenaAmbiguityResolution::solveAmbiguity(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
        measurementsPerTrack,
    Optional_cuts<track_container_t, traj_t, holder_t> optionalCuts) const {
  ACTS_INFO("Number of tracks before Ambiguty Resolution: " << tracks.size());
  std::vector<std::map<std::size_t, TrackFeatures>>
      trackFeaturesMaps;  // vector of trackFeaturesMaps. where each
                          // trackFeaturesMap contains the number of
                          // hits/hole/outliers for each detector in a track.
  std::vector<double> trackScore =
      simpleScore(tracks, trackFeaturesMaps, optionalCuts);

  std::vector<std::size_t> cleanTracks =
      getCleanedOutTracks(trackScore, trackFeaturesMaps, measurementsPerTrack);

  ACTS_VERBOSE("Number of clean tracks: " << cleanTracks.size());
  ACTS_VERBOSE("Min score: " << m_cfg.minScore);

  std::vector<int> goodTracks;
  std::size_t iTrack = 0;
  for (const auto& track : tracks) {
    if (std::find(cleanTracks.begin(), cleanTracks.end(), iTrack) !=
        cleanTracks.end()) {
      if (trackScore[iTrack] >= m_cfg.minScore) {
        goodTracks.push_back(track.index());
      }
    }
    iTrack++;
  }

  ACTS_INFO("Number of Good tracks: " << goodTracks.size());
  return goodTracks;
}

}  // namespace Acts
