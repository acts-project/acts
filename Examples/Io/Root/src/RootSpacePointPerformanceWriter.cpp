// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootSpacePointPerformanceWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/StripModulePairing.hpp"
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

template <typename Range1, typename Range2>
bool hasCommon(Range1&& a, Range2&& b) {
  for (const auto& elementA : a) {
    for (const auto& elementB : b) {
      if (elementA == elementB) {
        return true;
      }
    }
  }
  return false;
}

}  // namespace

RootSpacePointPerformanceWriter::RootSpacePointPerformanceWriter(
    const Config& config, Acts::Logging::Level level)
    : WriterT(config.inputSpacePoints, "RootSpacePointPerformanceWriter",
              level),
      m_cfg(config) {
  // Input container for space points is already checked by base constructor
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing sim hits input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing measurement-to-hits map input collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Missing measurement-to-particles map input collection");
  }
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);

  m_stripModulePairMap = pairStripModules(
      *m_cfg.trackingGeometry, m_cfg.stripGeometrySelection, logger());

  // Setup ROOT File
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }

  m_outputFile->cd();

  m_fakeVsZ.emplace("fake_vs_z", "Space point fake ratio vs z",
                    std::array{m_cfg.zAxis});
  m_fakeVsR.emplace("fake_vs_r", "Space point fake ratio vs r",
                    std::array{m_cfg.rAxis});
  m_fakeVsEta.emplace("fake_vs_eta", "Space point fake ratio vs eta",
                      std::array{m_cfg.etaAxis});
  m_fakeVsPhi.emplace("fake_vs_phi", "Space point fake ratio vs phi",
                      std::array{m_cfg.phiAxis});

  m_effVsZ.emplace("efficiency_vs_z", "Space point efficiency vs z",
                   std::array{m_cfg.zAxis});
  m_effVsR.emplace("efficiency_vs_r", "Space point efficiency vs r",
                   std::array{m_cfg.rAxis});
  m_effVsEta.emplace("efficiency_vs_eta", "Space point efficiency vs eta",
                     std::array{m_cfg.etaAxis});
  m_effVsPhi.emplace("efficiency_vs_phi", "Space point efficiency vs phi",
                     std::array{m_cfg.phiAxis});
}

RootSpacePointPerformanceWriter::~RootSpacePointPerformanceWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootSpacePointPerformanceWriter::finalize() {
  // Close the file if it's yours
  m_outputFile->cd();

  ActsPlugins::toRoot(*m_fakeVsZ)->Write();
  ActsPlugins::toRoot(*m_fakeVsR)->Write();
  ActsPlugins::toRoot(*m_fakeVsEta)->Write();
  ActsPlugins::toRoot(*m_fakeVsPhi)->Write();

  ActsPlugins::toRoot(*m_effVsZ)->Write();
  ActsPlugins::toRoot(*m_effVsR)->Write();
  ActsPlugins::toRoot(*m_effVsEta)->Write();
  ActsPlugins::toRoot(*m_effVsPhi)->Write();

  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

ProcessCode RootSpacePointPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const SpacePointContainer& spacePoints) {
  const auto& particles = m_inputParticles(ctx);
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& measurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);
  const auto& measurementParticlesMap = m_inputMeasurementParticlesMap(ctx);

  const auto transformSecond =
      std::views::transform([](auto& p) -> auto& { return p.second; });
  const auto filterParticles = std::views::filter(
      [&](const SimBarcode& particle) { return particles.contains(particle); });

  std::lock_guard<std::mutex> lock(m_writeMutex);

  std::map<std::pair<Index, Index>, SpacePointIndex> trueSpacePoints;

  for (const auto& spacePoint : spacePoints) {
    const Acts::Vector3 position(spacePoint.x(), spacePoint.y(),
                                 spacePoint.z());
    const double z = position.z();
    const double r = Acts::VectorHelpers::perp(position);
    const double phi = Acts::VectorHelpers::phi(position);
    const double eta = Acts::VectorHelpers::eta(position);

    const auto& sourceLinks = spacePoint.sourceLinks();

    bool isFake = false;

    if (sourceLinks.size() == 2) {
      const auto measurementIndex1 =
          sourceLinks[0].get<IndexSourceLink>().index();
      const auto measurementIndex2 =
          sourceLinks[1].get<IndexSourceLink>().index();

      auto contributingParticles1 =
          toRange(measurementParticlesMap.equal_range(measurementIndex1)) |
          transformSecond | filterParticles;
      auto contributingParticles2 =
          toRange(measurementParticlesMap.equal_range(measurementIndex2)) |
          transformSecond | filterParticles;

      isFake = !hasCommon(contributingParticles1, contributingParticles2);

      if (!isFake) {
        trueSpacePoints.emplace(
            std::pair<Index, Index>{measurementIndex1, measurementIndex2},
            spacePoint.index());
      }
    }

    m_fakeVsZ->fill({z}, isFake);
    m_fakeVsR->fill({r}, isFake);
    m_fakeVsEta->fill({eta}, isFake);
    m_fakeVsPhi->fill({phi}, isFake);
  }

  for (const auto [module1, module2] : m_stripModulePairMap) {
    const Acts::Surface* surface1 =
        m_cfg.trackingGeometry->findSurface(module1);
    const Acts::Surface* surface2 =
        m_cfg.trackingGeometry->findSurface(module2);

    if (surface1 == nullptr || surface2 == nullptr) {
      ACTS_WARNING("Could not find surfaces for modules " << module1 << " and "
                                                          << module2);
      continue;
    }

    const auto sourceLinks1 =
        measurements.orderedIndices().equal_range(module1);
    const auto sourceLinks2 =
        measurements.orderedIndices().equal_range(module2);

    for (const auto& sourceLink1 : toRange(sourceLinks1)) {
      const auto measurementIndex1 = sourceLink1.index();
      auto contributingParticles1 =
          toRange(measurementParticlesMap.equal_range(measurementIndex1)) |
          transformSecond | filterParticles;

      for (const auto& sourceLink2 : toRange(sourceLinks2)) {
        const auto measurementIndex2 = sourceLink2.index();
        auto contributingParticles2 =
            toRange(measurementParticlesMap.equal_range(measurementIndex2)) |
            transformSecond | filterParticles;

        if (!hasCommon(contributingParticles1, contributingParticles2)) {
          continue;
        }

        // Find the contributing simulated hits
        const auto hitRange1 =
            makeRange(measurementSimHitsMap.equal_range(measurementIndex1));
        // Use average truth in the case of multiple contributing sim hits
        const auto [local1, position1, dir1] = averageSimHits(
            ctx.geoContext, *surface1, simHits, hitRange1, logger());

        const double z = position1.z();
        const double r = Acts::VectorHelpers::perp(position1);
        const double phi = Acts::VectorHelpers::phi(position1);
        const double eta = Acts::VectorHelpers::eta(position1.head<3>());

        const bool isReconstructed =
            trueSpacePoints.contains(std::pair<Index, Index>{
                measurementIndex1, measurementIndex2}) ||
            trueSpacePoints.contains(
                std::pair<Index, Index>{measurementIndex2, measurementIndex1});

        m_effVsZ->fill({z}, isReconstructed);
        m_effVsR->fill({r}, isReconstructed);
        m_effVsEta->fill({eta}, isReconstructed);
        m_effVsPhi->fill({phi}, isReconstructed);
      }
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
