// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Performance/VertexPerformanceWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MultiIndex.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ios>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <stdexcept>
#include <utility>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {
struct AlgorithmContext;
}  // namespace ActsExamples

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
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Collection with particles missing");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing input track particles matching");
  }
  if (m_cfg.inputParticleTrackMatching.empty()) {
    throw std::invalid_argument("Missing input particle track matching");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);

  if (m_cfg.useTracks) {
    m_inputTracks.initialize(m_cfg.inputTracks);
  }

  // Setup ROOT I/O
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + path + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  } else {
    // I/O parameters.
    m_outputTree->Branch("event_nr", &m_eventNr);

    // Branches related to the 4D vertex position
    m_outputTree->Branch("truthX", &m_truthX);
    m_outputTree->Branch("truthY", &m_truthY);
    m_outputTree->Branch("truthZ", &m_truthZ);
    m_outputTree->Branch("truthT", &m_truthT);

    m_outputTree->Branch("recoX", &m_recoX);
    m_outputTree->Branch("recoY", &m_recoY);
    m_outputTree->Branch("recoZ", &m_recoZ);
    m_outputTree->Branch("recoT", &m_recoT);

    m_outputTree->Branch("resX", &m_resX);
    m_outputTree->Branch("resY", &m_resY);
    m_outputTree->Branch("resZ", &m_resZ);
    m_outputTree->Branch("resT", &m_resT);

    m_outputTree->Branch("pullX", &m_pullX);
    m_outputTree->Branch("pullY", &m_pullY);
    m_outputTree->Branch("pullZ", &m_pullZ);
    m_outputTree->Branch("pullT", &m_pullT);

    m_outputTree->Branch("covXX", &m_covXX);
    m_outputTree->Branch("covYY", &m_covYY);
    m_outputTree->Branch("covZZ", &m_covZZ);
    m_outputTree->Branch("covTT", &m_covTT);
    m_outputTree->Branch("covXY", &m_covXY);
    m_outputTree->Branch("covXZ", &m_covXZ);
    m_outputTree->Branch("covXT", &m_covXT);
    m_outputTree->Branch("covYZ", &m_covYZ);
    m_outputTree->Branch("covYT", &m_covYT);
    m_outputTree->Branch("covZT", &m_covZT);

    m_outputTree->Branch("sumPt2", &m_sumPt2);

    // Branches related to track momenta at vertex
    m_outputTree->Branch("trk_truthPhi", &m_truthPhi);
    m_outputTree->Branch("trk_truthTheta", &m_truthTheta);
    m_outputTree->Branch("trk_truthQOverP", &m_truthQOverP);

    m_outputTree->Branch("trk_recoPhi", &m_recoPhi);
    m_outputTree->Branch("trk_recoPhiFitted", &m_recoPhiFitted);
    m_outputTree->Branch("trk_recoTheta", &m_recoTheta);
    m_outputTree->Branch("trk_recoThetaFitted", &m_recoThetaFitted);
    m_outputTree->Branch("trk_recoQOverP", &m_recoQOverP);
    m_outputTree->Branch("trk_recoQOverPFitted", &m_recoQOverPFitted);

    m_outputTree->Branch("trk_resPhi", &m_resPhi);
    m_outputTree->Branch("trk_resPhiFitted", &m_resPhiFitted);
    m_outputTree->Branch("trk_resTheta", &m_resTheta);
    m_outputTree->Branch("trk_resThetaFitted", &m_resThetaFitted);
    m_outputTree->Branch("trk_resQOverP", &m_resQOverP);
    m_outputTree->Branch("trk_resQOverPFitted", &m_resQOverPFitted);
    m_outputTree->Branch("trk_momOverlap", &m_momOverlap);
    m_outputTree->Branch("trk_momOverlapFitted", &m_momOverlapFitted);

    m_outputTree->Branch("trk_pullPhi", &m_pullPhi);
    m_outputTree->Branch("trk_pullPhiFitted", &m_pullPhiFitted);
    m_outputTree->Branch("trk_pullTheta", &m_pullTheta);
    m_outputTree->Branch("trk_pullThetaFitted", &m_pullThetaFitted);
    m_outputTree->Branch("trk_pullQOverP", &m_pullQOverP);
    m_outputTree->Branch("trk_pullQOverPFitted", &m_pullQOverPFitted);

    m_outputTree->Branch("trk_weight", &m_trkWeight);

    m_outputTree->Branch("nTracksTruthVtx", &m_nTracksOnTruthVertex);
    m_outputTree->Branch("nTracksRecoVtx", &m_nTracksOnRecoVertex);

    m_outputTree->Branch("trkVtxMatch", &m_trackVtxMatchFraction);

    m_outputTree->Branch("nTrueVtx", &m_nTrueVtx);
    m_outputTree->Branch("nRecoVtx", &m_nRecoVtx);
    m_outputTree->Branch("nVtxReconstructable", &m_nVtxReconstructable);
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
    const AlgorithmContext& ctx, const std::vector<Acts::Vertex>& vertices) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_nRecoVtx = vertices.size();

  ACTS_DEBUG("Number of reco vertices in event: " << m_nRecoVtx);

  // Read truth input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);
  // TODO
  // const auto& particleTrackMatching = m_inputParticleTrackMatching(ctx);

  // Get number of generated true primary vertices
  m_nTrueVtx = getNumberOfTruePriVertices(particles);

  ACTS_VERBOSE("Number of truth particles in event : " << particles.size());
  ACTS_VERBOSE("Number of truth primary vertices : " << m_nTrueVtx);

  const ConstTrackContainer* tracks = nullptr;
  SimParticleContainer associatedParticles;

  // Get the event number
  m_eventNr = ctx.eventNumber;

  auto findParticle = [&](const auto& track) -> std::optional<SimParticle> {
    // Get the truth-matched particle
    auto imatched = trackParticleMatching.find(track.index());
    if (imatched == trackParticleMatching.end() ||
        !imatched->second.particle.has_value()) {
      ACTS_DEBUG("No truth particle associated with this track, index = "
                 << track.index() << " tip index = " << track.tipIndex());
      return {};
    }

    const auto& particleMatch = imatched->second;

    auto iparticle = particles.find(particleMatch.particle->value());
    if (iparticle == particles.end()) {
      ACTS_DEBUG(
          "Truth particle found but not monitored with this track, index = "
          << track.index() << " tip index = " << track.tipIndex()
          << " and this barcode = " << particleMatch.particle->value());
      return {};
    }

    return *iparticle;
  };

  if (m_cfg.useTracks) {
    tracks = &m_inputTracks(ctx);

    for (const auto& track : *tracks) {
      if (!track.hasReferenceSurface()) {
        ACTS_DEBUG("No reference surface on this track, index = "
                   << track.index() << " tip index = " << track.tipIndex());
        continue;
      }

      if (auto particle = findParticle(track); particle.has_value()) {
        associatedParticles.insert(*particle);
      }
    }
  } else {
    // if not using tracks, then all truth particles are associated with the
    // vertex
    associatedParticles = particles;
  }

  // Get number of track-associated true primary vertices
  m_nVtxReconstructable =
      getNumberOfReconstructableVertices(associatedParticles);

  ACTS_INFO("Number of reconstructed tracks : "
            << ((tracks != nullptr) ? tracks->size() : 0));
  ACTS_INFO("Number of reco track-associated truth particles in event : "
            << associatedParticles.size());
  ACTS_INFO("Maximum number of reconstructible primary vertices : "
            << m_nVtxReconstructable);

  // We compare the reconstructed momenta to the true momenta at the vertex. For
  // this, we propagate the reconstructed tracks to the PCA of the true vertex
  // position. Setting up propagator:
  Acts::EigenStepper<> stepper(m_cfg.bField);
  using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
  auto propagator = std::make_shared<Propagator>(stepper);
  // Setting the geometry/magnetic field context for the event
  Acts::PropagatorOptions pOptions(ctx.geoContext, ctx.magFieldContext);

  // Loop over reconstructed vertices and see if they can be matched to a true
  // vertex.
  for (const auto& vtx : vertices) {
    // Reconstructed tracks that contribute to the reconstructed vertex
    const auto& tracksAtVtx = vtx.tracks();
    // Input tracks matched to `tracksAtVtx`
    std::vector<std::uint32_t> trackIndices;

    // Containers for storing truth particles and truth vertices that contribute
    // to the reconstructed vertex
    SimParticleContainer particleAtVtx;
    SimBarcodeContainer contributingTruthVertices;

    if (m_cfg.useTracks) {
      for (const auto& trk : tracksAtVtx) {
        // Track parameters before the vertex fit
        const Acts::BoundTrackParameters& origTrack =
            *trk.originalParams.as<Acts::BoundTrackParameters>();

        bool foundMatchingParticle = false;

        // Finding the matching parameters in the container of all track
        // parameters. This allows us to identify the corresponding particle.
        // TODO this should not be necessary if the tracks at vertex would keep
        // this information
        for (const auto& inputTrk : *tracks) {
          const auto& params = inputTrk.parameters();

          if (origTrack.parameters() == params) {
            trackIndices.push_back(inputTrk.index());

            auto particleOpt = findParticle(inputTrk);
            if (!particleOpt.has_value()) {
              continue;
            }
            const auto& particle = *particleOpt;

            particleAtVtx.insert(particle);
            SimBarcode vtxId = SimBarcode(particle.particleId())
                                   .setParticle(0)
                                   .setGeneration(0)
                                   .setSubParticle(0);
            contributingTruthVertices.insert(vtxId);

            foundMatchingParticle = true;

            break;
          }
        }

        if (!foundMatchingParticle) {
          ACTS_VERBOSE("Track has no matching truth particle.");
        }
      }

      if (tracksAtVtx.size() != trackIndices.size()) {
        ACTS_ERROR("Not all tracks at vertex have a matching input track.");
      }
    } else {
      for (const auto& particle : particles) {
        SimBarcode vtxId = SimBarcode(particle.particleId())
                               .setParticle(0)
                               .setGeneration(0)
                               .setSubParticle(0);
        contributingTruthVertices.insert(vtxId);
      }
    }

    // Find true vertex that contributes most to the reconstructed vertex
    std::map<SimBarcode, std::uint32_t> fmap;
    for (SimBarcode vtxId : contributingTruthVertices) {
      fmap[vtxId]++;
    }
    std::uint32_t maxOccurrence = 0;
    SimBarcode maxOccurrenceId;
    for (const auto& [vtxId, count] : fmap) {
      if (count > maxOccurrence) {
        maxOccurrenceId = vtxId;
        maxOccurrence = count;
      }
    }

    // Count number of reconstructible tracks on truth vertex
    std::uint32_t nTracksOnTruthVertex = 0;
    for (const auto& particle : associatedParticles) {
      SimBarcode vtxId = SimBarcode(particle.particleId())
                             .setParticle(0)
                             .setGeneration(0)
                             .setSubParticle(0);
      if (vtxId == maxOccurrenceId) {
        ++nTracksOnTruthVertex;
      }
    }

    // Get number of contributing tracks (i.e., tracks with a weight above
    // threshold)
    auto weightHighEnough = [this](const auto& trkAtVtx) {
      return trkAtVtx.trackWeight > m_cfg.minTrkWeight;
    };
    unsigned int nTracksOnRecoVertex =
        std::count_if(tracksAtVtx.begin(), tracksAtVtx.end(), weightHighEnough);
    // Match reconstructed and truth vertex if the tracks of the truth vertex
    // make up at least minTrackVtxMatchFraction of the tracks at the
    // reconstructed vertex.
    double trackVtxMatchFraction =
        (m_cfg.useTracks ? (double)fmap[maxOccurrenceId] / nTracksOnRecoVertex
                         : 1.0);
    if (trackVtxMatchFraction > m_cfg.minTrackVtxMatchFraction) {
      // Get references to inner vectors where all track variables corresponding
      // to the current vertex will be saved
      auto& innerTruthPhi = m_truthPhi.emplace_back();
      auto& innerTruthTheta = m_truthTheta.emplace_back();
      auto& innerTruthQOverP = m_truthQOverP.emplace_back();

      auto& innerRecoPhi = m_recoPhi.emplace_back();
      auto& innerRecoTheta = m_recoTheta.emplace_back();
      auto& innerRecoQOverP = m_recoQOverP.emplace_back();

      auto& innerRecoPhiFitted = m_recoPhiFitted.emplace_back();
      auto& innerRecoThetaFitted = m_recoThetaFitted.emplace_back();
      auto& innerRecoQOverPFitted = m_recoQOverPFitted.emplace_back();

      auto& innerResPhi = m_resPhi.emplace_back();
      auto& innerResTheta = m_resTheta.emplace_back();
      auto& innerResQOverP = m_resQOverP.emplace_back();

      auto& innerResPhiFitted = m_resPhiFitted.emplace_back();
      auto& innerResThetaFitted = m_resThetaFitted.emplace_back();
      auto& innerResQOverPFitted = m_resQOverPFitted.emplace_back();

      auto& innerMomOverlap = m_momOverlap.emplace_back();
      auto& innerMomOverlapFitted = m_momOverlapFitted.emplace_back();

      auto& innerPullPhi = m_pullPhi.emplace_back();
      auto& innerPullTheta = m_pullTheta.emplace_back();
      auto& innerPullQOverP = m_pullQOverP.emplace_back();

      auto& innerPullPhiFitted = m_pullPhiFitted.emplace_back();
      auto& innerPullThetaFitted = m_pullThetaFitted.emplace_back();
      auto& innerPullQOverPFitted = m_pullQOverPFitted.emplace_back();

      auto& innerTrkWeight = m_trkWeight.emplace_back();

      // Helper function for computing the pull
      auto pull =
          [&](const Acts::ActsScalar& diff, const Acts::ActsScalar& variance,
              const std::string& variableStr, const bool& afterFit = true) {
            if (variance <= 0) {
              std::string tempStr;
              if (afterFit) {
                tempStr = "after";
              } else {
                tempStr = "before";
              }
              ACTS_WARNING("Nonpositive variance "
                           << tempStr << " vertex fit: Var(" << variableStr
                           << ") = " << variance << " <= 0.");
              return std::numeric_limits<double>::quiet_NaN();
            }
            double std = std::sqrt(variance);
            return diff / std;
          };

      const auto& truthVertexParticle =
          *std::find_if(associatedParticles.begin(), associatedParticles.end(),
                        [&](const auto& particle) {
                          return particle.particleId() == maxOccurrenceId;
                        });

      const Acts::ActsVector<4>& truePos = truthVertexParticle.fourPosition();

      // Write vertex truth based information
      {
        m_truthX.push_back(truePos[Acts::FreeIndices::eFreePos0]);
        m_truthY.push_back(truePos[Acts::FreeIndices::eFreePos1]);
        m_truthZ.push_back(truePos[Acts::FreeIndices::eFreePos2]);
        m_truthT.push_back(truePos[Acts::FreeIndices::eFreeTime]);

        m_recoX.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos0]);
        m_recoY.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos1]);
        m_recoZ.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos2]);
        m_recoT.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreeTime]);

        const Acts::ActsVector<4> diffPos = vtx.fullPosition() - truePos;
        m_resX.push_back(diffPos[Acts::FreeIndices::eFreePos0]);
        m_resY.push_back(diffPos[Acts::FreeIndices::eFreePos1]);
        m_resZ.push_back(diffPos[Acts::FreeIndices::eFreePos2]);
        m_resT.push_back(diffPos[Acts::FreeIndices::eFreeTime]);

        Acts::ActsScalar varX = vtx.fullCovariance()(
            Acts::FreeIndices::eFreePos0, Acts::FreeIndices::eFreePos0);
        Acts::ActsScalar varY = vtx.fullCovariance()(
            Acts::FreeIndices::eFreePos1, Acts::FreeIndices::eFreePos1);
        Acts::ActsScalar varZ = vtx.fullCovariance()(
            Acts::FreeIndices::eFreePos2, Acts::FreeIndices::eFreePos2);
        Acts::ActsScalar varTime = vtx.fullCovariance()(
            Acts::FreeIndices::eFreeTime, Acts::FreeIndices::eFreeTime);
        m_pullX.push_back(
            pull(diffPos[Acts::FreeIndices::eFreePos0], varX, "X"));
        m_pullY.push_back(
            pull(diffPos[Acts::FreeIndices::eFreePos1], varY, "Y"));
        m_pullZ.push_back(
            pull(diffPos[Acts::FreeIndices::eFreePos2], varZ, "Z"));
        m_pullT.push_back(
            pull(diffPos[Acts::FreeIndices::eFreeTime], varTime, "T"));

        m_covXX.push_back(varX);
        m_covYY.push_back(varY);
        m_covZZ.push_back(varZ);
        m_covTT.push_back(varTime);
        m_covXY.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos0,
                                               Acts::FreeIndices::eFreePos1));
        m_covXZ.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos0,
                                               Acts::FreeIndices::eFreePos2));
        m_covXT.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos0,
                                               Acts::FreeIndices::eFreeTime));
        m_covYZ.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos1,
                                               Acts::FreeIndices::eFreePos2));
        m_covYT.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos1,
                                               Acts::FreeIndices::eFreeTime));
        m_covZT.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos2,
                                               Acts::FreeIndices::eFreeTime));

        double sumPt2 = 0;
        for (const auto& trk : tracksAtVtx) {
          if (trk.trackWeight > m_cfg.minTrkWeight) {
            double pt = trk.originalParams.as<Acts::BoundTrackParameters>()
                            ->transverseMomentum();
            sumPt2 += pt * pt;
          }
        }
        m_sumPt2.push_back(sumPt2);

        m_nTracksOnTruthVertex.push_back(nTracksOnTruthVertex);
        m_nTracksOnRecoVertex.push_back(nTracksOnRecoVertex);

        m_trackVtxMatchFraction.push_back(trackVtxMatchFraction);
      }

      // Write vertex track based information
      {
        // Saving the reconstructed/truth momenta. The reconstructed momenta
        // are taken at the PCA to the truth vertex position -> we need to
        // perform a propagation.

        // Perigee at the true vertex position
        const std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
            Acts::Surface::makeShared<Acts::PerigeeSurface>(truePos.head(3));
        // Lambda for propagating the tracks to the PCA
        auto propagateToVtx = [&](const Acts::BoundTrackParameters& params)
            -> std::optional<Acts::BoundTrackParameters> {
          auto intersection =
              perigeeSurface
                  ->intersect(ctx.geoContext, params.position(ctx.geoContext),
                              params.direction(), Acts::BoundaryCheck(false))
                  .closest();
          pOptions.direction = Acts::Direction::fromScalarZeroAsPositive(
              intersection.pathLength());

          auto result =
              propagator->propagate(params, *perigeeSurface, pOptions);
          if (!result.ok()) {
            ACTS_ERROR("Propagation to true vertex position failed.");
            return std::nullopt;
          }
          auto& paramsAtVtx = *result->endParameters;
          return std::make_optional(paramsAtVtx);
        };

        for (const auto& [trkAtVtx, trkIndex] :
             Acts::zip(tracksAtVtx, trackIndices)) {
          const auto trk = tracks->getTrack(trkIndex);

          innerTrkWeight.push_back(trkAtVtx.trackWeight);

          auto particleOpt = findParticle(trk);
          if (!particleOpt.has_value()) {
            continue;
          }
          const auto& particle = *particleOpt;

          const auto& trueUnitDir = particle.direction();
          Acts::ActsVector<3> trueMom;
          trueMom.head(2) = Acts::makePhiThetaFromDirection(trueUnitDir);
          trueMom[2] = particle.qOverP();
          innerTruthPhi.push_back(trueMom[0]);
          innerTruthTheta.push_back(trueMom[1]);
          innerTruthQOverP.push_back(trueMom[2]);

          // Save track parameters before the vertex fit
          const auto paramsAtVtx = propagateToVtx(
              *trkAtVtx.originalParams.as<Acts::BoundTrackParameters>());
          if (paramsAtVtx != std::nullopt) {
            Acts::ActsVector<3> recoMom =
                paramsAtVtx->parameters().segment(Acts::eBoundPhi, 3);
            const Acts::ActsMatrix<3, 3>& momCov =
                paramsAtVtx->covariance()->block<3, 3>(Acts::eBoundPhi,
                                                       Acts::eBoundPhi);
            innerRecoPhi.push_back(recoMom[0]);
            innerRecoTheta.push_back(recoMom[1]);
            innerRecoQOverP.push_back(recoMom[2]);

            Acts::ActsVector<3> diffMom = recoMom - trueMom;
            // Accounting for the periodicity of phi. We overwrite the
            // previously computed value for better readability.
            diffMom[0] = Acts::detail::difference_periodic(
                recoMom(0), trueMom(0), 2 * M_PI);
            innerResPhi.push_back(diffMom[0]);
            innerResTheta.push_back(diffMom[1]);
            innerResQOverP.push_back(diffMom[2]);

            innerPullPhi.push_back(
                pull(diffMom[0], momCov(0, 0), "phi", false));
            innerPullTheta.push_back(
                pull(diffMom[1], momCov(1, 1), "theta", false));
            innerPullQOverP.push_back(
                pull(diffMom[2], momCov(2, 2), "q/p", false));

            const auto& recoUnitDir = paramsAtVtx->direction();
            double overlap = trueUnitDir.dot(recoUnitDir);
            innerMomOverlap.push_back(overlap);
          }

          // Save track parameters after the vertex fit
          const auto paramsAtVtxFitted = propagateToVtx(trkAtVtx.fittedParams);
          if (paramsAtVtxFitted != std::nullopt &&
              trkAtVtx.trackWeight > m_cfg.minTrkWeight) {
            Acts::ActsVector<3> recoMomFitted =
                paramsAtVtxFitted->parameters().segment(Acts::eBoundPhi, 3);
            const Acts::ActsMatrix<3, 3>& momCovFitted =
                paramsAtVtxFitted->covariance()->block<3, 3>(Acts::eBoundPhi,
                                                             Acts::eBoundPhi);
            innerRecoPhiFitted.push_back(recoMomFitted[0]);
            innerRecoThetaFitted.push_back(recoMomFitted[1]);
            innerRecoQOverPFitted.push_back(recoMomFitted[2]);

            Acts::ActsVector<3> diffMomFitted = recoMomFitted - trueMom;
            // Accounting for the periodicity of phi. We overwrite the
            // previously computed value for better readability.
            diffMomFitted[0] = Acts::detail::difference_periodic(
                recoMomFitted(0), trueMom(0), 2 * M_PI);
            innerResPhiFitted.push_back(diffMomFitted[0]);
            innerResThetaFitted.push_back(diffMomFitted[1]);
            innerResQOverPFitted.push_back(diffMomFitted[2]);

            innerPullPhiFitted.push_back(
                pull(diffMomFitted[0], momCovFitted(0, 0), "phi"));
            innerPullThetaFitted.push_back(
                pull(diffMomFitted[1], momCovFitted(1, 1), "theta"));
            innerPullQOverPFitted.push_back(
                pull(diffMomFitted[2], momCovFitted(2, 2), "q/p"));

            const auto& recoUnitDirFitted = paramsAtVtxFitted->direction();
            double overlapFitted = trueUnitDir.dot(recoUnitDirFitted);
            innerMomOverlapFitted.push_back(overlapFitted);
          }
        }
      }
    }
  }  // end loop vertices

  // fill the variables
  m_outputTree->Fill();

  m_truthX.clear();
  m_truthY.clear();
  m_truthZ.clear();
  m_truthT.clear();
  m_recoX.clear();
  m_recoY.clear();
  m_recoZ.clear();
  m_recoT.clear();
  m_resX.clear();
  m_resY.clear();
  m_resZ.clear();
  m_resT.clear();
  m_pullX.clear();
  m_pullY.clear();
  m_pullZ.clear();
  m_pullT.clear();
  m_covXX.clear();
  m_covYY.clear();
  m_covZZ.clear();
  m_covTT.clear();
  m_covXY.clear();
  m_covXZ.clear();
  m_covXT.clear();
  m_covYZ.clear();
  m_covYT.clear();
  m_covZT.clear();
  m_sumPt2.clear();
  m_truthPhi.clear();
  m_truthTheta.clear();
  m_truthQOverP.clear();
  m_recoPhi.clear();
  m_recoPhiFitted.clear();
  m_recoTheta.clear();
  m_recoThetaFitted.clear();
  m_recoQOverP.clear();
  m_recoQOverPFitted.clear();
  m_resPhi.clear();
  m_resPhiFitted.clear();
  m_resTheta.clear();
  m_resThetaFitted.clear();
  m_resQOverP.clear();
  m_resQOverPFitted.clear();
  m_momOverlap.clear();
  m_momOverlapFitted.clear();
  m_pullPhi.clear();
  m_pullPhiFitted.clear();
  m_pullTheta.clear();
  m_pullThetaFitted.clear();
  m_pullQOverP.clear();
  m_pullQOverPFitted.clear();

  m_trkWeight.clear();

  m_nTracksOnTruthVertex.clear();
  m_nTracksOnRecoVertex.clear();

  m_trackVtxMatchFraction.clear();

  return ProcessCode::SUCCESS;
}
