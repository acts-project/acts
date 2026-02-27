// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMeasurementPerformanceWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <ios>
#include <ranges>
#include <stdexcept>

#include <TEfficiency.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

namespace ActsExamples {

namespace {

template <typename Begin, typename End>
auto toRange(const std::pair<Begin, End>& pair) {
  return std::ranges::subrange(pair.first, pair.second);
}

}  // namespace

RootMeasurementPerformanceWriter::RootMeasurementPerformanceWriter(
    const Config& config, Acts::Logging::Level level)
    : WriterT(config.inputMeasurements, "RootMeasurementPerformanceWriter",
              level),
      m_cfg(config) {
  // Input container for measurements is already checked by base constructor
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing measurement-to-hits map input collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Missing measurement-to-particles map input collection");
  }
  if (m_cfg.inputSimHitMeasurementsMap.empty()) {
    throw std::invalid_argument(
        "Missing hits-to-measurements map input collection");
  }

  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_inputSimHitMeasurementsMap.initialize(m_cfg.inputSimHitMeasurementsMap);

  // Setup ROOT File
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }

  m_outputFile->cd();

  m_nContributingSimHits.emplace("nContributingSimHits",
                                 "Number of contributing simulated hits",
                                 std::array{m_cfg.countAxis});
  m_nContributingParticles.emplace("nContributingParticles",
                                   "Number of contributing particles",
                                   std::array{m_cfg.countAxis});
  m_purity.emplace("purity", "Measurements purity",
                   std::array{m_cfg.purityAxis});

  m_effVsZ.emplace("efficiency_vs_z", "Measurement efficiency vs z",
                   std::array{m_cfg.zAxis});
  m_effVsR.emplace("efficiency_vs_r", "Measurement efficiency vs r",
                   std::array{m_cfg.rAxis});
  m_effVsEta.emplace("efficiency_vs_eta", "Measurement efficiency vs eta",
                     std::array{m_cfg.etaAxis});
  m_effVsPhi.emplace("efficiency_vs_phi", "Measurement efficiency vs phi",
                     std::array{m_cfg.phiAxis});
}

RootMeasurementPerformanceWriter::~RootMeasurementPerformanceWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootMeasurementPerformanceWriter::finalize() {
  // Close the file if it's yours
  m_outputFile->cd();

  ActsPlugins::toRoot(*m_nContributingSimHits)->Write();
  ActsPlugins::toRoot(*m_nContributingParticles)->Write();
  ActsPlugins::toRoot(*m_purity)->Write();

  ActsPlugins::toRoot(*m_effVsZ)->Write();
  ActsPlugins::toRoot(*m_effVsR)->Write();
  ActsPlugins::toRoot(*m_effVsEta)->Write();
  ActsPlugins::toRoot(*m_effVsPhi)->Write();

  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

ProcessCode RootMeasurementPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const MeasurementContainer& measurements) {
  const auto& simHits = m_inputSimHits(ctx);
  const auto& measurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);
  const auto& measurementParticlesMap = m_inputMeasurementParticlesMap(ctx);
  const auto& simHitMeasurementsMap = m_inputSimHitMeasurementsMap(ctx);

  std::lock_guard<std::mutex> lock(m_writeMutex);

  for (const auto& measurement : measurements) {
    const auto contributingSimHits =
        toRange(measurementSimHitsMap.equal_range(measurement.index()));
    const auto contributingParticles =
        toRange(measurementParticlesMap.equal_range(measurement.index()));

    const double nContributingSimHits = contributingSimHits.size();
    const double nContributingParticles = contributingParticles.size();
    const double purity = 1.0 / nContributingParticles;

    m_nContributingSimHits->fill({nContributingSimHits});
    m_nContributingParticles->fill({nContributingParticles});
    m_purity->fill({purity});
  }

  for (const auto& simHit : simHits) {
    const auto relatedMeasurements =
        toRange(simHitMeasurementsMap.equal_range(simHit.index()));
    const std::size_t nRelatedMeasurements = relatedMeasurements.size();
    const bool efficient = nRelatedMeasurements > 0;

    const Acts::Vector3& position = simHit.position();
    const Acts::Vector3 normalizedPosition = simHit.position().normalized();
    const double z = position.z();
    const double r = Acts::VectorHelpers::perp(position);
    const double phi = Acts::VectorHelpers::phi(normalizedPosition);
    const double eta = Acts::VectorHelpers::eta(normalizedPosition);

    m_effVsZ->fill({z}, efficient);
    m_effVsR->fill({r}, efficient);
    m_effVsEta->fill({eta}, efficient);
    m_effVsPhi->fill({phi}, efficient);
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
