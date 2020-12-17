// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "SeedingPerformanceWriter.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <stdexcept>
#include <unordered_map>

#include <TFile.h>

using HitParticlesMap = ActsExamples::IndexMultimap<ActsFatras::Barcode>;

ActsExamples::SeedingPerformanceWriter::SeedingPerformanceWriter(
    ActsExamples::SeedingPerformanceWriter::Config cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputSeeds, "SeedingPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl) {
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
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
  auto path = joinPaths(m_cfg.outputDir, m_cfg.outputFilename);
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
  ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);
  const auto hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);

  for (auto& regionVec : seedVector) {
    nSeeds += regionVec.size();
    for (size_t i = 0; i < regionVec.size(); i++) {
      const Acts::Seed<SimSpacePoint>* seed = &regionVec[i];
      ProtoTrack ptrack{seed->sp()[0]->measurementIndex(),
                        seed->sp()[1]->measurementIndex(),
                        seed->sp()[2]->measurementIndex()};
      std::vector<ParticleHitCount> particleHitCounts;
      identifyContributingParticles(hitParticlesMap, ptrack, particleHitCounts);

      if (particleHitCounts.size() == 1) {
        auto prt = particleHitCounts.at(0);
        if (prt.hitCount == 3) {
          auto it = truthCount.try_emplace(prt.particleId, 0u).first;
          it->second += 1;
          nMatchedSeeds++;
        }
      }
    }
  }

  int nMatchedParticles = 0;
  int nDuplicatedParticles = 0;
  // Fill the effeciency and fake rate plots
  for (const auto& particle : particles) {
    const auto it1 = truthCount.find(particle.particleId());
    bool isMatched = false;
    int nMatchedSeedsForParticle = 0;
    if (it1 != truthCount.end()) {
      isMatched = true;
      nMatchedParticles++;
      nMatchedSeedsForParticle = it1->second;
      if (nMatchedSeedsForParticle > 1) {
        nDuplicatedParticles++;
      }
    }
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
