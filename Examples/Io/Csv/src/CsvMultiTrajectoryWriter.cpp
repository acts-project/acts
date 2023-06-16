// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvMultiTrajectoryWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/MultiIndex.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace ActsExamples {
class IndexSourceLink;
}  // namespace ActsExamples

using namespace ActsExamples;

CsvMultiTrajectoryWriter::CsvMultiTrajectoryWriter(
    const CsvMultiTrajectoryWriter::Config& config, Acts::Logging::Level level)
    : WriterT<TrajectoriesContainer>(config.inputTrajectories,
                                     "CsvMultiTrajectoryWriter", level),
      m_cfg(config) {
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectories collection");
  }

  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
}

ProcessCode CsvMultiTrajectoryWriter::writeT(
    const AlgorithmContext& context,
    const TrajectoriesContainer& trajectories) {
  // open per-event file
  std::string path =
      perEventFilepath(m_cfg.outputDir, m_cfg.fileName, context.eventNumber);
  std::ofstream mos(path, std::ofstream::out | std::ofstream::trunc);
  if (!mos) {
    throw std::ios_base::failure("Could not open '" + path + "' to write");
  }

  const auto& hitParticlesMap = m_inputMeasurementParticlesMap(context);

  std::unordered_map<Acts::MultiTrajectoryTraits::IndexType, trackInfo> infoMap;

  // Counter of truth-matched reco tracks
  using RecoTrackInfo = std::pair<trackInfo, size_t>;
  std::map<ActsFatras::Barcode, std::vector<RecoTrackInfo>> matched;

  size_t trackId = 0;
  for (const auto& traj : trajectories) {
    // The trajectory entry indices and the multiTrajectory
    const auto& trackTips = traj.tips();
    const auto& mj = traj.multiTrajectory();
    if (trackTips.empty()) {
      ACTS_WARNING("Empty multiTrajectory.");
      continue;
    }

    // Loop over all trajectories in a multiTrajectory
    for (auto trackTip : trackTips) {
      // Collect the trajectory summary info
      auto trajState =
          Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
      // Reco track selection
      //@TODO: add interface for applying others cuts on reco tracks:
      // -> pT, d0, z0, detector-specific hits/holes number cut
      if (trajState.nMeasurements < m_cfg.nMeasurementsMin) {
        continue;
      }

      // Check if the reco track has fitted track parameters
      if (not traj.hasTrackParameters(trackTip)) {
        ACTS_WARNING(
            "No fitted track parameters for trajectory with entry index = "
            << trackTip);
        continue;
      }

      // Get the majority truth particle to this track
      std::vector<ParticleHitCount> particleHitCount;
      identifyContributingParticles(hitParticlesMap, traj, trackTip,
                                    particleHitCount);
      if (m_cfg.onlyTruthMatched && particleHitCount.empty()) {
        ACTS_WARNING(
            "No truth particle associated with this trajectory with entry "
            "index = "
            << trackTip);
        continue;
      }

      // Requirement on the pT of the track
      const auto& momentum = traj.trackParameters(trackTip).momentum();
      const auto pT = Acts::VectorHelpers::perp(momentum);
      if (pT < m_cfg.ptMin) {
        continue;
      }
      size_t nMajorityHits = 0;
      ActsFatras::Barcode majorityParticleId;
      if (!particleHitCount.empty()) {
        // Get the majority particle counts
        majorityParticleId = particleHitCount.front().particleId;
        // n Majority hits
        nMajorityHits = particleHitCount.front().hitCount;
      }
      // track info
      trackInfo toAdd;
      toAdd.trackId = trackId;
      toAdd.particleId = majorityParticleId;
      toAdd.nStates = trajState.nStates;
      toAdd.nMajorityHits = nMajorityHits;
      toAdd.nMeasurements = trajState.nMeasurements;
      toAdd.nOutliers = trajState.nOutliers;
      toAdd.nHoles = trajState.nHoles;
      toAdd.nSharedHits = trajState.nSharedHits;
      toAdd.chi2Sum = trajState.chi2Sum;
      toAdd.NDF = trajState.NDF;
      toAdd.truthMatchProb = toAdd.nMajorityHits * 1. / trajState.nMeasurements;
      toAdd.fittedParameters = &traj.trackParameters(trackTip);
      toAdd.trackType = "unknown";

      mj.visitBackwards(trackTip, [&](const auto& state) {
        if (state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag) ==
            true) {
          auto sl =
              state.getUncalibratedSourceLink().template get<IndexSourceLink>();
          auto hitIndex = sl.index();
          toAdd.measurementsID.insert(toAdd.measurementsID.begin(), hitIndex);
        }
      });

      // Check if the trajectory is matched with truth.
      if (toAdd.truthMatchProb >= m_cfg.truthMatchProbMin) {
        matched[toAdd.particleId].push_back({toAdd, toAdd.trackId});
      } else {
        toAdd.trackType = "fake";
      }

      infoMap[toAdd.trackId] = toAdd;

      trackId++;
    }  // end of one trajectory
  }    // end of multi-trajectory

  // Find duplicates
  std::unordered_set<size_t> listGoodTracks;
  for (auto& [particleId, matchedTracks] : matched) {
    std::sort(matchedTracks.begin(), matchedTracks.end(),
              [](const RecoTrackInfo& lhs, const RecoTrackInfo& rhs) {
                // sort by nMajorityHits
                if (lhs.first.nMajorityHits != rhs.first.nMajorityHits) {
                  return (lhs.first.nMajorityHits > rhs.first.nMajorityHits);
                }
                // sort by nOutliers
                if (lhs.first.nOutliers != rhs.first.nOutliers) {
                  return (lhs.first.nOutliers < rhs.first.nOutliers);
                }
                // sort by chi2
                return (lhs.first.chi2Sum < rhs.first.chi2Sum);
              });

    listGoodTracks.insert(matchedTracks.front().first.trackId);
  }

  // write csv header
  mos << "track_id,particleId,"
      << "nStates,nMajorityHits,nMeasurements,nOutliers,nHoles,nSharedHits,"
      << "chi2,ndf,chi2/ndf,"
      << "pT,eta,phi,"
      << "truthMatchProbability,"
      << "good/duplicate/fake,"
      << "Hits_ID";

  mos << '\n';
  mos << std::setprecision(m_cfg.outputPrecision);

  // good/duplicate/fake = 0/1/2
  for (auto& [id, trajState] : infoMap) {
    if (listGoodTracks.find(id) != listGoodTracks.end()) {
      trajState.trackType = "good";
    } else if (trajState.trackType != "fake") {
      trajState.trackType = "duplicate";
    }

    // write the track info
    mos << trajState.trackId << ",";
    mos << trajState.particleId << ",";
    mos << trajState.nStates << ",";
    mos << trajState.nMajorityHits << ",";
    mos << trajState.nMeasurements << ",";
    mos << trajState.nOutliers << ",";
    mos << trajState.nHoles << ",";
    mos << trajState.nSharedHits << ",";
    mos << trajState.chi2Sum << ",";
    mos << trajState.NDF << ",";
    mos << trajState.chi2Sum * 1.0 / trajState.NDF << ",";
    mos << Acts::VectorHelpers::perp(trajState.fittedParameters->momentum())
        << ",";
    mos << Acts::VectorHelpers::eta(trajState.fittedParameters->momentum())
        << ",";
    mos << Acts::VectorHelpers::phi(trajState.fittedParameters->momentum())
        << ",";
    mos << trajState.truthMatchProb << ",";
    mos << trajState.trackType << ",";
    mos << "\"[";
    for (auto& ID : trajState.measurementsID) {
      mos << ID << ",";
    }
    mos << "]\"";
    mos << '\n';
  }

  return ProcessCode::SUCCESS;
}
