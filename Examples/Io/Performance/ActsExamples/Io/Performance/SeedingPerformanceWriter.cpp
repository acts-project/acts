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
#include <set>
#include <stdexcept>
#include <TFile.h>

using HitParticlesMap = ActsExamples::IndexMultimap<ActsFatras::Barcode>;

ActsExamples::SeedingPerformanceWriter::SeedingPerformanceWriter(
    ActsExamples::SeedingPerformanceWriter::Config cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputSeeds, "SeedingPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl),
      m_nTotalSeeds(0),
      m_nTotalMatchedSeeds(0),
      m_nTotalParticles(0),
      m_nTotalMatchedParticles(0),
      m_nTotalDuplicatedParticles(0) {
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing input seeds collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  if (m_cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.inputHitParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = m_cfg.outputFilename;
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (not m_outputFile) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }
  // initialize the plot tools
  m_effPlotTool.book(m_effPlotCache);
}

ActsExamples::SeedingPerformanceWriter::~SeedingPerformanceWriter() {
  m_effPlotTool.clear(m_effPlotCache);

  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::SeedingPerformanceWriter::endRun() {
  float eff = float(m_nTotalMatchedParticles) / m_nTotalParticles;
  float fakeRate = float(m_nTotalSeeds - m_nTotalMatchedSeeds) / m_nTotalSeeds;
  float duplicationRate =
      float(m_nTotalDuplicatedParticles) / m_nTotalMatchedParticles;

  ACTS_INFO("Efficiency (nMatchedParticles / nAllParticles) = " << eff);
  ACTS_INFO("Fake rate (nUnMatchedSeeds / nAllSeeds) =" << fakeRate);
  ACTS_INFO(
      "Duplication rate (nDuplicatedMatchedParticles / nMatchedParticles) ="
      << duplicationRate);

  if (m_outputFile) {
    m_outputFile->cd();
    m_effPlotTool.write(m_effPlotCache);
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

std::set<ActsFatras::Barcode>
ActsExamples::SeedingPerformanceWriter::identifySharedParticles(
    const HitParticlesMap& hitParticlesMap,
    const Acts::Seed<SimSpacePoint>* seed) const {
  auto particles0 = getTruthParticles(hitParticlesMap, seed->sp()[0]->Id());
  auto particles1 = getTruthParticles(hitParticlesMap, seed->sp()[1]->Id());
  auto particles2 = getTruthParticles(hitParticlesMap, seed->sp()[2]->Id());

  std::set<ActsFatras::Barcode> tmp;
  set_intersection(particles0.begin(), particles0.end(), particles1.begin(),
                   particles1.end(), std::inserter(tmp, tmp.end()));

  std::set<ActsFatras::Barcode> prtsInCommon;
  set_intersection(particles2.begin(), particles2.end(), tmp.begin(), tmp.end(),
                   std::inserter(prtsInCommon, prtsInCommon.end()));
  return prtsInCommon;
}

// get truth particles that are a part of this space point
std::vector<ActsFatras::Barcode>
ActsExamples::SeedingPerformanceWriter::getTruthParticles(
    const HitParticlesMap& hitParticlesMap, const std::size_t hit_id) const {
  std::vector<ActsExamples::ParticleHitCount> particleHitCount;
  std::vector<ActsFatras::Barcode> particles;
  for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hit_id))) {
    auto particleId = hitParticle.second;
    // search for existing particle in the existing hit counts
    auto isSameParticle = [=](const ParticleHitCount& phc) {
      return (phc.particleId == particleId);
    };
    auto it = std::find_if(particleHitCount.begin(), particleHitCount.end(),
                           isSameParticle);
    if (it == particleHitCount.end()) {
      particles.push_back(particleId);
    }
  }
  return particles;
}

ActsExamples::ProcessCode ActsExamples::SeedingPerformanceWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<std::vector<Acts::Seed<SimSpacePoint>>>& seedVector) {
  // Read truth particles from input collection
  const auto& particles =
      ctx.eventStore.get<ActsExamples::SimParticleContainer>(
          m_cfg.inputParticles);

  size_t nSeeds = 0;
  size_t nMatchedSeeds = 0;
  // Map from particles to how many times they were successfully found by a seed
  std::unordered_map<ActsFatras::Barcode, std::size_t> truthCount;
  const HitParticlesMap hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputHitParticlesMap);

  for (auto& regionVec : seedVector) {
    nSeeds += regionVec.size();
    for (size_t i = 0; i < regionVec.size(); i++) {
      const Acts::Seed<SimSpacePoint>* seed = &regionVec[i];

      auto prtsInCommon = identifySharedParticles(hitParticlesMap, seed);
      if (prtsInCommon.size() > 0) {
        for (const auto& prt : prtsInCommon) {
          auto it = truthCount.try_emplace(prt, 0u).first;
          it->second += 1;
        }
        nMatchedSeeds++;
      }
    }
  }

  int nMatchedParticles = 0;
  int nDuplicatedParticles = 0;
  // Fill the effeciency and fake rate plots
  for (const auto& particle : particles) {
    const auto it1 = truthCount.find(particle.particleId());
    bool isMatched = false;
    bool isDuplicated = false;
    int nMatchedSeedsForParticle = 0;
    if (it1 != truthCount.end()) {
      isMatched = true;
      nMatchedParticles++;
      nMatchedSeedsForParticle = truthCount[particle.particleId()];
      if (nMatchedSeedsForParticle > 1)
        isDuplicated = true;
    }
    if (isDuplicated)
      nDuplicatedParticles++;
    m_effPlotTool.fill(m_effPlotCache, particle, isMatched);
  }

  ACTS_INFO("Number of seeds: " << nSeeds);
  m_nTotalSeeds += nSeeds;
  m_nTotalMatchedSeeds += nMatchedSeeds;
  m_nTotalParticles += particles.size();
  m_nTotalMatchedParticles += nMatchedParticles;
  m_nTotalDuplicatedParticles += nDuplicatedParticles;

  return ProcessCode::SUCCESS;
}
