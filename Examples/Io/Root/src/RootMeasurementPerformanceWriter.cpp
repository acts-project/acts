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

  m_measurementContributingHits.emplace(
      "measurement_n_contributing_hits",
      "Measurement number of contributing simulated hits",
      std::array{m_cfg.countAxis});
  m_measurementContributingParticles.emplace(
      "measurement_n_contributing_particles",
      "Measurement number of contributing particles",
      std::array{m_cfg.countAxis});
  m_measurementPurity.emplace("measurement_purity", "Measurements purity",
                              std::array{m_cfg.purityAxis});
  m_measurementClassification.emplace(
      "measurement_classification", "Measurement classification",
      std::array{Acts::Experimental::AxisVariant(
          Acts::Experimental::BoostRegularAxis{4, 0, 4, "Classification"})});

  m_hitEffVsZ.emplace("hit_efficiency_vs_z",
                      "Hit to measurement efficiency vs z",
                      std::array{m_cfg.zAxis});
  m_hitEffVsR.emplace("hit_efficiency_vs_r",
                      "Hit to measurement efficiency vs r",
                      std::array{m_cfg.rAxis});
  m_hitEffVsEta.emplace("hit_efficiency_vs_eta",
                        "Hit to measurement efficiency vs eta",
                        std::array{m_cfg.etaAxis});
  m_hitEffVsPhi.emplace("hit_efficiency_vs_phi",
                        "Hit to measurement efficiency vs phi",
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

  ActsPlugins::toRoot(*m_measurementContributingHits)->Write();
  ActsPlugins::toRoot(*m_measurementContributingParticles)->Write();
  ActsPlugins::toRoot(*m_measurementPurity)->Write();
  ActsPlugins::toRoot(*m_measurementClassification)->Write();

  ActsPlugins::toRoot(*m_hitEffVsZ)->Write();
  ActsPlugins::toRoot(*m_hitEffVsR)->Write();
  ActsPlugins::toRoot(*m_hitEffVsEta)->Write();
  ActsPlugins::toRoot(*m_hitEffVsPhi)->Write();

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
    const MeasurementClassification classification = [&]() {
      if (nContributingParticles == 0) {
        return MeasurementClassification::Fake;
      } else if (purity > m_cfg.matchingRatio) {
        return MeasurementClassification::Matched;
      } else {
        return MeasurementClassification::Merged;
      }
    }();

    m_measurementContributingHits->fill({nContributingSimHits});
    m_measurementContributingParticles->fill({nContributingParticles});
    m_measurementPurity->fill({purity});
    m_measurementClassification->fill({static_cast<double>(classification)});
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

    m_hitEffVsZ->fill({z}, efficient);
    m_hitEffVsR->fill({r}, efficient);
    m_hitEffVsEta->fill({eta}, efficient);
    m_hitEffVsPhi->fill({phi}, efficient);
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
