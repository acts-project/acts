// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <cstddef>
#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <RtypesCore.h>
#include <TFile.h>
#include <TTree.h>

struct ActsExamples::TrackFinderPerformanceWriter::Impl {
  Config cfg;

  ReadDataHandle<SimParticleContainer> inputParticles;
  ReadDataHandle<HitParticlesMap> inputMeasurementParticlesMap;
  ReadDataHandle<TrackParticleMatching> inputTrackParticleMatching;

  TFile* file = nullptr;

  // per-track tree
  TTree* trkTree = nullptr;
  std::mutex trkMutex;
  // track identification
  ULong64_t trkEventId = 0;
  ULong64_t trkTrackId = 0;
  // track content
  // number of hits on track
  UShort_t trkNumHits = 0;
  // number of particles contained in the track
  UShort_t trkNumParticles = 0;
  // track particle content; for each contributing particle, largest first
  std::vector<ULong64_t> trkParticleId;
  // total number of hits generated by this particle
  std::vector<UShort_t> trkParticleNumHitsTotal;
  // number of hits within this track
  std::vector<UShort_t> trkParticleNumHitsOnTrack;

  // per-particle tree
  TTree* prtTree = nullptr;
  std::mutex prtMutex;
  // particle identification
  ULong64_t prtEventId = 0;
  ULong64_t prtParticleId = 0;
  Int_t prtParticleType = 0;
  // particle kinematics
  // vertex position in mm
  float prtVx = 0, prtVy = 0, prtVz = 0;
  // vertex time in ns
  float prtVt = 0;
  // particle momentum at production in GeV
  float prtPx = 0, prtPy = 0, prtPz = 0;
  // particle mass in GeV
  float prtM = 0;
  // particle charge in e
  float prtQ = 0;
  // particle reconstruction
  UShort_t prtNumHits = 0;  // number of hits for this particle
  UShort_t prtNumTracks =
      0;  // number of tracks this particle was reconstructed in
  UShort_t prtNumTracksMajority =
      0;  // number of tracks reconstructed as majority
  // extra logger reference for the logging macros
  const Acts::Logger& _logger;

  Impl(TrackFinderPerformanceWriter* parent, Config&& c, const Acts::Logger& l)
      : cfg(std::move(c)),
        inputParticles{parent, "InputParticles"},
        inputMeasurementParticlesMap{parent, "InputMeasurementParticlesMap"},
        inputTrackParticleMatching{parent, "InputTrackParticleMatching"},
        _logger(l) {
    if (cfg.inputTracks.empty()) {
      throw std::invalid_argument("Missing track input collection");
    }
    if (cfg.inputParticles.empty()) {
      throw std::invalid_argument("Missing particles input collection");
    }
    if (cfg.inputMeasurementParticlesMap.empty()) {
      throw std::invalid_argument("Missing hit-particles map input collection");
    }
    if (cfg.inputTrackParticleMatching.empty()) {
      throw std::invalid_argument(
          "Missing proto track-particle matching input collection");
    }
    if (cfg.filePath.empty()) {
      throw std::invalid_argument("Missing output filename");
    }

    inputParticles.initialize(cfg.inputParticles);
    inputMeasurementParticlesMap.initialize(cfg.inputMeasurementParticlesMap);
    inputTrackParticleMatching.initialize(cfg.inputTrackParticleMatching);

    // the output file can not be given externally since TFile accesses to the
    // same file from multiple threads are unsafe.
    // must always be opened internally
    file = TFile::Open(cfg.filePath.c_str(), cfg.fileMode.c_str());
    if (file == nullptr) {
      throw std::invalid_argument("Could not open '" + cfg.filePath + "'");
    }

    // construct trees
    trkTree = new TTree(cfg.treeNameTracks.c_str(), cfg.treeNameTracks.c_str());
    trkTree->SetDirectory(file);
    trkTree->Branch("event_id", &trkEventId);
    trkTree->Branch("track_id", &trkTrackId);
    trkTree->Branch("size", &trkNumHits);
    trkTree->Branch("nparticles", &trkNumParticles);
    trkTree->Branch("particle_id", &trkParticleId);
    trkTree->Branch("particle_nhits_total", &trkParticleNumHitsTotal);
    trkTree->Branch("particle_nhits_on_track", &trkParticleNumHitsOnTrack);
    prtTree =
        new TTree(cfg.treeNameParticles.c_str(), cfg.treeNameParticles.c_str());
    prtTree->SetDirectory(file);
    prtTree->Branch("event_id", &prtEventId);
    prtTree->Branch("particle_id", &prtParticleId);
    prtTree->Branch("particle_type", &prtParticleType);
    prtTree->Branch("vx", &prtVx);
    prtTree->Branch("vy", &prtVy);
    prtTree->Branch("vz", &prtVz);
    prtTree->Branch("vt", &prtVt);
    prtTree->Branch("px", &prtPx);
    prtTree->Branch("py", &prtPy);
    prtTree->Branch("pz", &prtPz);
    prtTree->Branch("m", &prtM);
    prtTree->Branch("q", &prtQ);
    prtTree->Branch("nhits", &prtNumHits);
    prtTree->Branch("ntracks", &prtNumTracks);
    prtTree->Branch("ntracks_majority", &prtNumTracksMajority);
  }

  const Acts::Logger& logger() const { return _logger; }

  void write(std::uint64_t eventId, const ConstTrackContainer& tracks,
             const SimParticleContainer& particles,
             const HitParticlesMap& hitParticlesMap,
             const TrackParticleMatching& trackParticleMatching) {
    const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);

    // How often a particle was reconstructed.
    std::unordered_map<ActsFatras::Barcode, std::size_t> reconCount;
    reconCount.reserve(particles.size());
    // How often a particle was reconstructed as the majority particle.
    std::unordered_map<ActsFatras::Barcode, std::size_t> majorityCount;
    majorityCount.reserve(particles.size());

    // write per-track performance measures
    {
      std::lock_guard<std::mutex> guardTrk(trkMutex);
      for (auto track : tracks) {
        // Get the truth-matched particle
        auto imatched = trackParticleMatching.find(track.index());
        if (imatched == trackParticleMatching.end()) {
          ACTS_DEBUG(
              "No truth particle associated with this proto track, index = "
              << track.index());
          continue;
        }
        const auto& particleMatch = imatched->second;

        if (particleMatch.particle.has_value()) {
          SimBarcode majorityParticleId = particleMatch.particle.value();

          auto it = majorityCount.try_emplace(majorityParticleId, 0u).first;
          it->second += 1;

          // Find the truth particle via the barcode
          if (auto ip = particles.find(majorityParticleId);
              ip == particles.end()) {
            ACTS_WARNING(
                "Majority particle not found in the particles collection.");
          }
        }

        for (const auto& hc : particleMatch.contributingParticles) {
          auto it = reconCount.try_emplace(hc.particleId, 0u).first;
          it->second += 1;
        }

        trkEventId = eventId;
        trkTrackId = track.index();
        trkNumHits = track.nMeasurements();
        trkNumParticles = particleMatch.contributingParticles.size();
        trkParticleId.clear();
        trkParticleNumHitsTotal.clear();
        trkParticleNumHitsOnTrack.clear();
        for (const auto& phc : particleMatch.contributingParticles) {
          trkParticleId.push_back(phc.particleId.value());
          // count total number of hits for this particle
          auto trueParticleHits =
              makeRange(particleHitsMap.equal_range(phc.particleId.value()));
          trkParticleNumHitsTotal.push_back(trueParticleHits.size());
          trkParticleNumHitsOnTrack.push_back(phc.hitCount);
        }

        trkTree->Fill();
      }
    }

    // write per-particle performance measures
    {
      std::lock_guard<std::mutex> guardPrt(trkMutex);
      for (const auto& particle : particles) {
        // find all hits for this particle
        auto hits =
            makeRange(particleHitsMap.equal_range(particle.particleId()));

        // identification
        prtEventId = eventId;
        prtParticleId = particle.particleId().value();
        prtParticleType = particle.pdg();
        // kinematics
        prtVx = particle.position().x() / Acts::UnitConstants::mm;
        prtVy = particle.position().y() / Acts::UnitConstants::mm;
        prtVz = particle.position().z() / Acts::UnitConstants::mm;
        prtVt = particle.time() / Acts::UnitConstants::mm;
        const auto p = particle.absoluteMomentum() / Acts::UnitConstants::GeV;
        prtPx = p * particle.direction().x();
        prtPy = p * particle.direction().y();
        prtPz = p * particle.direction().z();
        prtM = particle.mass() / Acts::UnitConstants::GeV;
        prtQ = particle.charge() / Acts::UnitConstants::e;
        // reconstruction
        prtNumHits = hits.size();
        auto nt = reconCount.find(particle.particleId());
        prtNumTracks = (nt != reconCount.end()) ? nt->second : 0u;
        auto nm = majorityCount.find(particle.particleId());
        prtNumTracksMajority = (nm != majorityCount.end()) ? nm->second : 0u;

        prtTree->Fill();
      }
    }
  }

  /// Write everything to disk and close the file.
  void close() {
    if (file == nullptr) {
      ACTS_ERROR("Output file is not available");
      return;
    }
    file->Write();
    file->Close();
  }
};

ActsExamples::TrackFinderPerformanceWriter::TrackFinderPerformanceWriter(
    ActsExamples::TrackFinderPerformanceWriter::Config config,
    Acts::Logging::Level level)
    : WriterT(config.inputTracks, "TrackFinderPerformanceWriter", level),
      m_impl(std::make_unique<Impl>(this, std::move(config), logger())) {}

ActsExamples::TrackFinderPerformanceWriter::~TrackFinderPerformanceWriter() =
    default;

ActsExamples::ProcessCode ActsExamples::TrackFinderPerformanceWriter::writeT(
    const ActsExamples::AlgorithmContext& ctx,
    const ActsExamples::ConstTrackContainer& tracks) {
  const auto& particles = m_impl->inputParticles(ctx);
  const auto& hitParticlesMap = m_impl->inputMeasurementParticlesMap(ctx);
  const auto& trackParticleMatching = m_impl->inputTrackParticleMatching(ctx);
  m_impl->write(ctx.eventNumber, tracks, particles, hitParticlesMap,
                trackParticleMatching);
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode
ActsExamples::TrackFinderPerformanceWriter::finalize() {
  m_impl->close();
  return ProcessCode::SUCCESS;
}

const ActsExamples::TrackFinderPerformanceWriter::Config&
ActsExamples::TrackFinderPerformanceWriter::config() const {
  return m_impl->cfg;
}
