// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <stdexcept>

std::vector<bool> Acts::ScoreBasedAmbiguityResolution::getCleanedOutTracks(
    const std::vector<double>& trackScore,
    const std::vector<std::vector<TrackFeatures>>& trackFeaturesMaps,
    const std::vector<std::vector<measurementTuple>>& measurementsPerTrack)
    const {
  std::vector<bool> cleanTracks(measurementsPerTrack.size(), false);

  ACTS_VERBOSE("Cleaning tracks");

  if (trackScore.size() != measurementsPerTrack.size()) {
    throw std::invalid_argument(
        "Track score and measurementsPerTrack size mismatch");
  }

  std::size_t numberOfTracks = measurementsPerTrack.size();
  ACTS_DEBUG("Number of tracks: " << numberOfTracks);

  boost::container::flat_map<std::size_t,
                             boost::container::flat_set<std::size_t>>
      tracksPerMeasurement;

  // Removes bad tracks and counts computes the vector of tracks per
  // measurement.
  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    if (trackScore[iTrack] <= 0) {
      continue;
    }
    for (auto measurementTuples : measurementsPerTrack[iTrack]) {
      auto iMeasurement = std::get<0>(measurementTuples);
      tracksPerMeasurement[iMeasurement].insert(iTrack);
    }
  }

  enum TrackStateTypes {
    // A measurement not yet used in any other track
    UnsharedHit = 1,
    // A measurement shared with another track
    SharedHit = 2,
    // A hit that needs to be removed from the track
    RejectedHit = 3,
    // an outlier, to be copied in case
    Outlier = 4,
    // other trackstate types to be copied in case
    OtherTrackStateType = 5
  };

  std::vector<std::vector<std::size_t>> newMeasurements;
  // Loop over all tracks in the track container
  for (std::size_t iTrack = 0; iTrack < numberOfTracks; ++iTrack) {
    double track_score = trackScore[iTrack];
    ACTS_DEBUG("Track score: " << track_score);

    if (track_score <= 0) {
      ACTS_DEBUG("Track " << iTrack << " could not be accepted - low score");
      continue;
    }

    auto trackFeaturesMap = trackFeaturesMaps[iTrack];

    bool TrkCouldBeAccepted = true;

    // For tracks with shared hits, we need to check and remove bad hits

    std::vector<int> trackStateTypes(measurementsPerTrack[iTrack].size(),
                                     OtherTrackStateType);
    int index = 0;

    // Loop over all measurements of the track and for each hit a
    // trackStateTypes is assigned.
    for (auto measurementTuples : measurementsPerTrack[iTrack]) {
      auto iMeasurement = std::get<0>(measurementTuples);
      auto iVolume = std::get<1>(measurementTuples);
      auto isoutliner = std::get<2>(measurementTuples);

      auto volume_it = m_cfg.volumeMap.find(iVolume);

      if (volume_it == m_cfg.volumeMap.end()) {
        index++;
        continue;
      }

      auto detectorId = volume_it->second;
      auto detector_it = m_cfg.detectorMap.find(detectorId);
      auto detector = detector_it->second;

      if (isoutliner) {
        ACTS_VERBOSE("Measurement is outlier on a fitter track, copy it over");
        trackStateTypes[index] = Outlier;
        index++;
        continue;
      }
      if (tracksPerMeasurement[iMeasurement].size() == 1) {
        ACTS_VERBOSE("Measurement is not shared, copy it over");

        trackStateTypes[index] = UnsharedHit;

        index++;
        continue;
      }
      if (tracksPerMeasurement[iMeasurement].size() > 1) {
        ACTS_VERBOSE("Measurement is shared, copy it over");

        if (detector.sharedHitsFlag == true) {
          ACTS_VERBOSE("Measurement is shared, Reject it");
          trackStateTypes[index] = RejectedHit;
          index++;
          continue;
        }

        trackStateTypes[index] = SharedHit;

        index++;
        continue;
      }
    }

    std::vector<std::size_t> newMeasurementsPerTrack;
    std::size_t measurement = 0;
    std::size_t nshared = 0;

    // Loop over all measurements of the track and process them according to the
    // trackStateTypes and other conditions.
    // Good measurements are copied to the newMeasurementsPerTrack vector.
    for (std::size_t i = 0; i < trackStateTypes.size(); i++) {
      auto measurementTuples = measurementsPerTrack[iTrack][i];
      measurement = std::get<0>(measurementTuples);

      if (trackStateTypes[i] == RejectedHit) {
        ACTS_DEBUG("Dropping rejected hit");
      } else if (trackStateTypes[i] != SharedHit) {
        ACTS_DEBUG("Good TSOS, copy hit");
        newMeasurementsPerTrack.push_back(measurement);

        // a counter called nshared is used to keep track of the number of
        // shared hits accepted.
      } else if (nshared >= m_cfg.maxShared) {
        ACTS_DEBUG("Too many shared hit, drop it");
      }
      // If the track is shared, the hit is only accepted if the track has score
      // higher than the minimum score for shared tracks.
      else {
        ACTS_DEBUG("Try to recover shared hit ");
        if (tracksPerMeasurement[measurement].size() <
                m_cfg.maxSharedTracksPerMeasurement &&
            track_score > m_cfg.minScoreSharedTracks) {
          ACTS_DEBUG("Accepted hit shared with "
                     << tracksPerMeasurement[measurement].size() << " tracks");
          newMeasurementsPerTrack.push_back(measurement);
          nshared++;
        } else {
          ACTS_DEBUG("Rejected hit shared with "
                     << tracksPerMeasurement[measurement].size() << " tracks");
        }
      }
    }

    // Check if the track has enough hits to be accepted.
    if (newMeasurementsPerTrack.size() < 3) {
      TrkCouldBeAccepted = false;
      ACTS_DEBUG(std::endl
                 << "Track " << iTrack
                 << " could not be accepted - not enough hits");
      ACTS_DEBUG("Number of hits: " << measurementsPerTrack[iTrack].size());
      ACTS_DEBUG("Number of good hits: " << newMeasurementsPerTrack.size());
      continue;
    }

    // Check if the track has too many shared hits to be accepted.
    for (std::size_t detectorId = 0; detectorId < m_cfg.detectorMap.size();
         detectorId++) {
      auto detector_it = m_cfg.detectorMap.find(detectorId);
      auto detector = detector_it->second;
      if (trackFeaturesMap[detectorId].nSharedHits > detector.maxSharedHits) {
        TrkCouldBeAccepted = false;
        break;
      }
    }

    if (TrkCouldBeAccepted) {
      cleanTracks[iTrack] = true;
      newMeasurements.push_back(newMeasurementsPerTrack);
      ACTS_VERBOSE("Track " << iTrack << " is accepted");
      continue;
    }
  }

  return cleanTracks;
}
