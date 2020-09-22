// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "SeedingPerformanceWriter.hpp"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <numeric>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

namespace {
using Acts::VectorHelpers::eta;
using SimParticleContainer = ActsExamples::SimParticleContainer;
using ProtoTrackContainer = ActsExamples::ProtoTrackContainer;
using ProtoTrack = ActsExamples::ProtoTrack;
}  // namespace

ActsExamples::SeedingPerformanceWriter::SeedingPerformanceWriter(
    ActsExamples::SeedingPerformanceWriter::Config cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputSeeds, "SeedingPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl) {
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing input seeds collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  if (m_cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = m_cfg.outputFilename;
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (not m_outputFile) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }
  // initialize the efficiency plots tool
  m_effPlotTool.book(m_effPlotCache);
}

ActsExamples::SeedingPerformanceWriter::~SeedingPerformanceWriter() {
  m_effPlotTool.clear(m_effPlotCache);

  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::SeedingPerformanceWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_effPlotTool.write(m_effPlotCache);
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

std::set<ActsFatras::Barcode>
ActsExamples::SeedingPerformanceWriter::identifySharedParticles(
    const Acts::Seed<SimSpacePoint>* seed) const {
  const SimSpacePoint* sp0 = seed->sp()[0];
  const SimSpacePoint* sp1 = seed->sp()[1];
  const SimSpacePoint* sp2 = seed->sp()[2];
  std::set<ActsFatras::Barcode> particles0;
  std::set<ActsFatras::Barcode> particles1;
  std::set<ActsFatras::Barcode> particles2;
  for (size_t i = 0; i < sp0->particles.size(); i++) {
    particles0.insert(sp0->particles[i].particleId);  // insert particle barcode
  }
  for (size_t i = 0; i < sp1->particles.size(); i++) {
    particles1.insert(sp1->particles[i].particleId);
  }
  for (size_t i = 0; i < sp2->particles.size(); i++) {
    particles2.insert(sp2->particles[i].particleId);
  }
  std::set<ActsFatras::Barcode> tmp;
  set_intersection(particles0.begin(), particles0.end(), particles1.begin(),
                   particles1.end(), std::inserter(tmp, tmp.end()));
  std::set<ActsFatras::Barcode> prtsInCommon;
  set_intersection(particles2.begin(), particles2.end(), tmp.begin(), tmp.end(),
                   std::inserter(prtsInCommon, prtsInCommon.end()));
  return prtsInCommon;
}

ActsExamples::ProcessCode ActsExamples::SeedingPerformanceWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<std::vector<Acts::Seed<SimSpacePoint>>>& seedVector) {
  // Read truth particles from input collection
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);
  size_t nSeeds = 0;
  // Set that helps keep track of the number of duplicate seeds. i.e. if there
  // are three seeds for one particle, this is counted as two duplicate seeds.
  std::set<ActsFatras::Barcode> particlesFoundBySeeds;

  // Map from particles to how many times they were successfully found by a seed
  std::unordered_map<ActsFatras::Barcode, std::size_t> truthCount;

  for (auto& regionVec : seedVector) {
    nSeeds += regionVec.size();
    for (size_t i = 0; i < regionVec.size(); i++) {
      const Acts::Seed<SimSpacePoint>* seed = &regionVec[i];

      std::set<ActsFatras::Barcode> prtsInCommon =
          identifySharedParticles(seed);
      if (prtsInCommon.size() > 0) {
        for (const auto& prt : prtsInCommon) {
          auto it = truthCount.try_emplace(prt, 0u).first;
          it->second += 1;
        }
      }
    }
  }
  ACTS_INFO("Number of seeds: " << nSeeds);

  // Fill the effeciency and fake rate plots
  for (const auto& particle : particles) {
    const auto it1 = truthCount.find(particle.particleId());
    bool isMatched = false;

    if (it1 != truthCount.end()) {
      isMatched = true;
    }

    m_effPlotTool.fill(m_effPlotCache, particle, isMatched);
  }

  return ProcessCode::SUCCESS;
}
