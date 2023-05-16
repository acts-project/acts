// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Performance/VertexPerformanceWriter.hpp"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::VertexPerformanceWriter::VertexPerformanceWriter(
    const ActsExamples::VertexPerformanceWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputVertices, "VertexPerformanceWriter", level),
      m_cfg(config) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  if (m_cfg.inputAllTruthParticles.empty()) {
    throw std::invalid_argument("Collection with all truth particles missing");
  }
  if (m_cfg.inputSelectedTruthParticles.empty()) {
    throw std::invalid_argument(
        "Collection with selected truth particles missing");
  }

  m_inputAllTruthParticles.initialize(m_cfg.inputAllTruthParticles);
  m_inputSelectedTruthParticles.initialize(m_cfg.inputSelectedTruthParticles);

  if (!m_cfg.inputAssociatedTruthParticles.empty()) {
    m_inputAssociatedTruthParticles.initialize(
        m_cfg.inputAssociatedTruthParticles);
    if (!m_cfg.inputTrackParameters.empty()) {
      m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
    } else {
      m_inputTrajectories.initialize(m_cfg.inputTrajectories);
    }
  } else {
    m_inputMeasurementParticlesMap.initialize(
        m_cfg.inputMeasurementParticlesMap);
    m_inputTrajectories.initialize(m_cfg.inputTrajectories);
  }

  m_inputTime.maybeInitialize(m_cfg.inputTime);

  // Setup ROOT I/O
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + path);
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  } else {
    // I/O parameters
    m_outputTree->Branch("resX", &m_resX);
    m_outputTree->Branch("resY", &m_resY);
    m_outputTree->Branch("resZ", &m_resZ);

    m_outputTree->Branch("pullX", &m_pullX);
    m_outputTree->Branch("pullY", &m_pullY);
    m_outputTree->Branch("pullZ", &m_pullZ);

    m_outputTree->Branch("recoX", &m_recoX);
    m_outputTree->Branch("recoY", &m_recoY);
    m_outputTree->Branch("recoZ", &m_recoZ);

    m_outputTree->Branch("truthX", &m_truthX);
    m_outputTree->Branch("truthY", &m_truthY);
    m_outputTree->Branch("truthZ", &m_truthZ);

    m_outputTree->Branch("covXX", &m_covXX);
    m_outputTree->Branch("covYY", &m_covYY);
    m_outputTree->Branch("covZZ", &m_covZZ);
    m_outputTree->Branch("covXY", &m_covXY);
    m_outputTree->Branch("covXZ", &m_covXZ);
    m_outputTree->Branch("covYZ", &m_covYZ);

    m_outputTree->Branch("trkVtxMatch", &m_trackVtxMatchFraction);
    m_outputTree->Branch("nRecoVtx", &m_nrecoVtx);
    m_outputTree->Branch("nTrueVtx", &m_ntrueVtx);
    m_outputTree->Branch("nVtxDetectorAcceptance", &m_nVtxDetAcceptance);
    m_outputTree->Branch("nVtxReconstructable", &m_nVtxReconstructable);
    m_outputTree->Branch("timeMS", &m_timeMS);
  }
}

ActsExamples::VertexPerformanceWriter::~VertexPerformanceWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::VertexPerformanceWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

int ActsExamples::VertexPerformanceWriter::getNumberOfReconstructableVertices(
    const SimParticleContainer& collection) const {
  // map for finding frequency
  std::map<int, int> fmap;

  std::vector<int> reconstructableTruthVertices;

  // traverse the array for frequency
  for (const auto& p : collection) {
    int secVtxId = p.particleId().vertexSecondary();
    if (secVtxId != 0) {
      // truthparticle from secondary vtx
      continue;
    }
    int priVtxId = p.particleId().vertexPrimary();
    fmap[priVtxId]++;
  }

  // iterate over the map
  for (auto it : fmap) {
    // Require at least 2 tracks
    if (it.second > 1) {
      reconstructableTruthVertices.push_back(it.first);
    }
  }

  return reconstructableTruthVertices.size();
}

int ActsExamples::VertexPerformanceWriter::getNumberOfTruePriVertices(
    const SimParticleContainer& collection) const {
  // Vector to store indices of all primary vertices
  std::set<int> allPriVtxIds;
  for (const auto& p : collection) {
    int priVtxId = p.particleId().vertexPrimary();
    int secVtxId = p.particleId().vertexSecondary();
    if (secVtxId != 0) {
      // truthparticle from secondary vtx
      continue;
    }
    // Insert to set, removing duplicates
    allPriVtxIds.insert(priVtxId);
  }
  // Size of set corresponds to total number of primary vertices
  return allPriVtxIds.size();
}

ActsExamples::ProcessCode ActsExamples::VertexPerformanceWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_nrecoVtx = vertices.size();

  ACTS_DEBUG("Number of reco vertices in event: " << m_nrecoVtx);

  // Read truth particle input collection
  const auto& allTruthParticles = m_inputAllTruthParticles(ctx);
  // Get number of generated true primary vertices
  m_ntrueVtx = getNumberOfTruePriVertices(allTruthParticles);

  ACTS_VERBOSE("Total number of generated truth particles in event : "
               << allTruthParticles.size());
  ACTS_VERBOSE(
      "Total number of generated truth primary vertices : " << m_ntrueVtx);

  // Read selected truth particle input collection
  const auto& selectedTruthParticles = m_inputSelectedTruthParticles(ctx);
  // Get number of detector-accepted true primary vertices
  m_nVtxDetAcceptance = getNumberOfTruePriVertices(selectedTruthParticles);

  ACTS_VERBOSE("Total number of selected truth particles in event : "
               << selectedTruthParticles.size());
  ACTS_VERBOSE("Total number of detector-accepted truth primary vertices : "
               << m_nVtxDetAcceptance);

  std::vector<Acts::BoundTrackParameters> trackParameters;
  SimParticleContainer associatedTruthParticles;

  if (!m_cfg.inputAssociatedTruthParticles.empty()) {
    if (!m_cfg.inputTrackParameters.empty()) {
      trackParameters = m_inputTrackParameters(ctx);
    } else {
      const auto& inputTrajectories = m_inputTrajectories(ctx);

      for (const auto& trajectories : inputTrajectories) {
        for (auto tip : trajectories.tips()) {
          if (!trajectories.hasTrackParameters(tip)) {
            continue;
          }
          trackParameters.push_back(trajectories.trackParameters(tip));
        }
      }
    }

    // Read track-associated truth particle input collection
    associatedTruthParticles = m_inputAssociatedTruthParticles(ctx);

    /*****************  Start x,y,z resolution plots here *****************/
    // Matching tracks at vertex to fitted tracks that are in turn matched
    // to truth particles. Match reco and true vtx if >50% of tracks match

    auto mismatchMsg = [&](auto level, const auto& extra) {
      ACTS_LOG(level,
               "Number of fitted tracks and associated truth particles do not "
               "match. ("
                   << trackParameters.size()
                   << " != " << associatedTruthParticles.size()
                   << ") Not able to match fitted tracks at reconstructed "
                      "vertex to truth vertex."
                   << extra);
    };

    if (associatedTruthParticles.size() < trackParameters.size()) {
      mismatchMsg(Acts::Logging::ERROR, " Switch to hit based truth matching.");
    } else if (associatedTruthParticles.size() > trackParameters.size()) {
      mismatchMsg(Acts::Logging::INFO,
                  " This is likely due to track efficiency < 1");
    }
  } else {
    // get active tips
    const auto& inputTrajectories = m_inputTrajectories(ctx);

    std::vector<ParticleHitCount> particleHitCounts;

    const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);

    for (const auto& trajectories : inputTrajectories) {
      for (auto tip : trajectories.tips()) {
        if (!trajectories.hasTrackParameters(tip)) {
          continue;
        }

        trackParameters.push_back(trajectories.trackParameters(tip));

        identifyContributingParticles(hitParticlesMap, trajectories, tip,
                                      particleHitCounts);
        ActsFatras::Barcode majorityParticleId =
            particleHitCounts.front().particleId;
        size_t nMajorityHits = particleHitCounts.front().hitCount;

        auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(
            trajectories.multiTrajectory(), tip);

        if (nMajorityHits * 1. / trajState.nMeasurements <
            m_cfg.truthMatchProbMin) {
          continue;
        }

        auto it = std::find_if(allTruthParticles.begin(),
                               allTruthParticles.end(), [&](const auto& tp) {
                                 return tp.particleId() == majorityParticleId;
                               });

        if (it == allTruthParticles.end()) {
          continue;
        }

        const auto& majorityParticle = *it;
        associatedTruthParticles.emplace_hint(associatedTruthParticles.end(),
                                              majorityParticle);
      }
    }
  }

  ACTS_VERBOSE(
      "Total number of reconstructed tracks : " << trackParameters.size());

  // Get number of track-associated true primary vertices
  m_nVtxReconstructable =
      getNumberOfReconstructableVertices(associatedTruthParticles);

  ACTS_INFO("Total number of reco track-associated truth particles in event : "
            << associatedTruthParticles.size());
  ACTS_INFO("Total number of reco track-associated truth primary vertices : "
            << m_nVtxReconstructable);

  // Loop over all reco vertices and find associated truth particles
  std::vector<SimParticleContainer> truthParticlesAtVtxContainer;
  for (const auto& vtx : vertices) {
    const auto tracks = vtx.tracks();
    // Store all associated truth particles to current vtx
    SimParticleContainer particleAtVtx;

    std::vector<int> contributingTruthVertices;

    for (const auto& trk : tracks) {
      Acts::BoundTrackParameters origTrack = *(trk.originalParams);

      // Find associated truth particle now
      int idx = 0;
      for (const auto& particle : associatedTruthParticles) {
        if (origTrack.parameters() == trackParameters[idx].parameters()) {
          particleAtVtx.insert(particleAtVtx.end(), particle);

          int priVtxId = particle.particleId().vertexPrimary();
          contributingTruthVertices.push_back(priVtxId);
        }
        idx++;
      }
    }  // end loop tracks

    // Now find true vtx with most matching tracks at reco vtx
    // and check if it contributes more than 50 of all tracks
    std::map<int, int> fmap;
    for (int priVtxId : contributingTruthVertices) {
      fmap[priVtxId]++;
    }
    int maxOccurrence = -1;
    int maxOccurrenceId = -1;
    for (auto it : fmap) {
      if (it.second > maxOccurrence) {
        maxOccurrence = it.second;
        maxOccurrenceId = it.first;
      }
    }

    // Match reco to truth vertex if at least 50% of tracks match
    double trackVtxMatchFraction =
        (double)fmap[maxOccurrenceId] / tracks.size();
    if (trackVtxMatchFraction > m_cfg.minTrackVtxMatchFraction) {
      for (const auto& particle : associatedTruthParticles) {
        int priVtxId = particle.particleId().vertexPrimary();
        int secVtxId = particle.particleId().vertexSecondary();

        if (secVtxId != 0) {
          // truthparticle from secondary vtx
          continue;
        }

        if (priVtxId == maxOccurrenceId) {
          // Vertex found, fill varibles
          const auto& truePos = particle.position();

          m_resX.push_back(vtx.position()[0] - truePos[0]);
          m_resY.push_back(vtx.position()[1] - truePos[1]);
          m_resZ.push_back(vtx.position()[2] - truePos[2]);

          auto pull = [&](int i) {
            double var = vtx.covariance()(i, i);
            if (var < 0) {
              ACTS_WARNING("var(" << i << ") = " << var << " < 0");
              return std::numeric_limits<double>::quiet_NaN();
            }
            double std = std::sqrt(var);
            if (std == 0) {
              ACTS_WARNING("std(" << i << ") = 0");
              return std::numeric_limits<double>::quiet_NaN();
            }
            return (vtx.position()[i] - truePos[i]) / std;
          };

          m_pullX.push_back(pull(0));
          m_pullY.push_back(pull(1));
          m_pullZ.push_back(pull(2));

          m_truthX.push_back(truePos[0]);
          m_truthY.push_back(truePos[1]);
          m_truthZ.push_back(truePos[2]);

          m_recoX.push_back(vtx.position()[0]);
          m_recoY.push_back(vtx.position()[1]);
          m_recoZ.push_back(vtx.position()[2]);

          m_covXX.push_back(vtx.covariance()(0, 0));
          m_covYY.push_back(vtx.covariance()(1, 1));
          m_covZZ.push_back(vtx.covariance()(2, 2));
          m_covXY.push_back(vtx.covariance()(0, 1));
          m_covXZ.push_back(vtx.covariance()(0, 2));
          m_covYZ.push_back(vtx.covariance()(1, 2));

          m_trackVtxMatchFraction.push_back(trackVtxMatchFraction);
          // Next vertex now
          break;
        }
      }
    }
  }  // end loop vertices

  // Retrieve and set reconstruction time
  if (!m_cfg.inputTime.empty()) {
    const auto& reconstructionTimeMS = m_inputTime(ctx);
    m_timeMS = reconstructionTimeMS;
  } else {
    m_timeMS = -1;
  }

  // fill the variables
  m_outputTree->Fill();

  m_resX.clear();
  m_resY.clear();
  m_resZ.clear();
  m_pullX.clear();
  m_pullY.clear();
  m_pullZ.clear();
  m_truthX.clear();
  m_truthY.clear();
  m_truthZ.clear();
  m_recoX.clear();
  m_recoY.clear();
  m_recoZ.clear();
  m_covXX.clear();
  m_covYY.clear();
  m_covZZ.clear();
  m_covXY.clear();
  m_covXZ.clear();
  m_covYZ.clear();
  m_trackVtxMatchFraction.clear();

  return ProcessCode::SUCCESS;
}
