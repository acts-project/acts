// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootSpacePointWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {

namespace {

Acts::Vector4 closestPointOnLine(const Acts::Vector4& a, const Acts::Vector4& b,
                                 const Acts::Vector3& p) {
  const Acts::Vector4 ab = b - a;

  const double t =
      (p - a.head<3>()).dot(ab.head<3>()) / ab.head<3>().squaredNorm();

  return a + t * ab;
}

}  // namespace

RootSpacePointWriter::RootSpacePointWriter(
    const RootSpacePointWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputSpacePoints, "RootSpacePointWriter", level),
      m_cfg(config) {
  // inputParticles is already checked by base constructor
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  if (m_cfg.inputSimHits.empty() !=
      m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Input sim hits and measurement-to-particles map must both be given or "
        "both be empty");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty() !=
      m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Input measurement-to-particles map and measurement-to-sim-hits map "
        "must both be given or both be empty");
  }

  m_inputSimHits.maybeInitialize(m_cfg.inputSimHits);
  m_inputMeasurementParticlesMap.maybeInitialize(
      m_cfg.inputMeasurementParticlesMap);
  m_inputMeasurementSimHitsMap.maybeInitialize(
      m_cfg.inputMeasurementSimHitsMap);

  // open root file and create the tree
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  // setup the branches
  m_outputTree->Branch("event_id", &m_eventId);
  m_outputTree->Branch("measurement_id", &m_measurementId1, "measurement_id/l");
  m_outputTree->Branch("geometry_id", &m_geometryId1, "geometry_id/l");
  m_outputTree->Branch("measurement_id_2", &m_measurementId2,
                       "measurement_id_2/l");
  m_outputTree->Branch("geometry_id_2", &m_geometryId2, "geometry_id_2/l");
  m_outputTree->Branch("x", &m_x);
  m_outputTree->Branch("y", &m_y);
  m_outputTree->Branch("z", &m_z);
  m_outputTree->Branch("t", &m_t);
  m_outputTree->Branch("r", &m_r);
  m_outputTree->Branch("var_r", &m_var_r);
  m_outputTree->Branch("var_z", &m_var_z);

  if (m_inputSimHits.isInitialized()) {
    m_outputTree->Branch("fake", &m_fake);

    m_outputTree->Branch("true_x", &m_trueX);
    m_outputTree->Branch("true_y", &m_trueY);
    m_outputTree->Branch("true_z", &m_trueZ);
    m_outputTree->Branch("true_t", &m_trueT);
    m_outputTree->Branch("true_r", &m_trueR);

    m_outputTree->Branch("residual_rphi", &m_residualRPhi);
    m_outputTree->Branch("residual_z", &m_residualZ);
    m_outputTree->Branch("residual_dca", &m_residualDCA);
  }
}

RootSpacePointWriter::~RootSpacePointWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootSpacePointWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_VERBOSE("Wrote hits to tree '" << m_cfg.treeName << "' in '"
                                      << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode RootSpacePointWriter::writeT(
    const AlgorithmContext& ctx, const SpacePointContainer& spacePoints) {
  // ensure exclusive access to tree/file while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  SimHitContainer emptySimHits{};
  MeasurementParticlesMap emptyMeasurementParticlesMap{};
  MeasurementSimHitsMap emptyMeasurementSimHitsMap{};

  const auto& simHits =
      m_inputSimHits.isInitialized() ? m_inputSimHits(ctx) : emptySimHits;
  const auto& measurementParticlesMap =
      m_inputMeasurementParticlesMap.isInitialized()
          ? m_inputMeasurementParticlesMap(ctx)
          : emptyMeasurementParticlesMap;
  const auto& measurementSimHitsMap =
      m_inputMeasurementSimHitsMap.isInitialized()
          ? m_inputMeasurementSimHitsMap(ctx)
          : emptyMeasurementSimHitsMap;

  // Get the event number
  m_eventId = ctx.eventNumber;

  for (const auto& sp : spacePoints) {
    const auto& sl1 = sp.sourceLinks()[0].get<IndexSourceLink>();
    m_measurementId1 = sl1.index();
    m_geometryId1 = sl1.geometryId().value();
    m_volumeId1 = static_cast<std::uint32_t>(sl1.geometryId().volume());
    m_layerId1 = static_cast<std::uint32_t>(sl1.geometryId().layer());
    m_surfaceId1 = static_cast<std::uint32_t>(sl1.geometryId().sensitive());
    m_extraId1 = static_cast<std::uint32_t>(sl1.geometryId().extra());

    if (sp.sourceLinks().size() == 2) {
      const auto& sl2 = sp.sourceLinks()[1].get<IndexSourceLink>();
      m_measurementId2 = sl2.index();
      m_geometryId2 = sl2.geometryId().value();
      m_volumeId2 = static_cast<std::uint32_t>(sl2.geometryId().volume());
      m_layerId2 = static_cast<std::uint32_t>(sl2.geometryId().layer());
      m_surfaceId2 = static_cast<std::uint32_t>(sl2.geometryId().sensitive());
      m_extraId2 = static_cast<std::uint32_t>(sl2.geometryId().extra());
    } else {
      m_measurementId2 = 0;
      m_geometryId2 = 0;
      m_volumeId2 = 0;
      m_layerId2 = 0;
      m_surfaceId2 = 0;
      m_extraId2 = 0;
    }

    // write sp position
    m_x = sp.x() / Acts::UnitConstants::mm;
    m_y = sp.y() / Acts::UnitConstants::mm;
    m_z = sp.z() / Acts::UnitConstants::mm;
    m_t = sp.time() / Acts::UnitConstants::mm;
    m_r = sp.r() / Acts::UnitConstants::mm;
    m_var_r = sp.varianceR() / Acts::square(Acts::UnitConstants::mm);
    m_var_z = sp.varianceZ() / Acts::square(Acts::UnitConstants::mm);

    if (m_inputSimHits.isInitialized()) {
      const Acts::Vector4 recoPosition(sp.x(), sp.y(), sp.z(), sp.time());

      const Acts::Surface* surface1 =
          m_cfg.trackingGeometry->findSurface(sl1.geometryId());
      const auto hitRange1 =
          makeRange(measurementSimHitsMap.equal_range(m_measurementId1));
      const auto [trueLocal1, truePosition1, trueDir1] = averageSimHits(
          ctx.geoContext, *surface1, simHits, hitRange1, logger());

      if (sp.sourceLinks().size() == 1) {
        m_fake = false;

        m_trueX = truePosition1.x() / Acts::UnitConstants::mm;
        m_trueY = truePosition1.y() / Acts::UnitConstants::mm;
        m_trueZ = truePosition1.z() / Acts::UnitConstants::mm;
        m_trueT = truePosition1.w() / Acts::UnitConstants::mm;
        m_trueR =
            Acts::VectorHelpers::perp(truePosition1) / Acts::UnitConstants::mm;

        m_residualRPhi = m_trueR * Acts::detail::difference_periodic(
                                       Acts::VectorHelpers::phi(recoPosition),
                                       Acts::VectorHelpers::phi(truePosition1),
                                       2 * std::numbers::pi);
        m_residualZ = m_trueZ - recoPosition.z() / Acts::UnitConstants::mm;
        m_residualDCA = Acts::fastHypot(m_residualRPhi, m_residualZ);
      } else {
        const auto [p1b, p1e] =
            measurementParticlesMap.equal_range(m_measurementId1);
        const auto [p2b, p2e] =
            measurementParticlesMap.equal_range(m_measurementId2);

        m_fake = true;
        for (auto it1 = p1b; it1 != p1e; ++it1) {
          for (auto it2 = p2b; it2 != p2e; ++it2) {
            if (it1->second == it2->second) {
              m_fake = false;
            }
          }
        }

        if (m_fake) {
          m_trueX = std::numeric_limits<float>::quiet_NaN();
          m_trueY = std::numeric_limits<float>::quiet_NaN();
          m_trueZ = std::numeric_limits<float>::quiet_NaN();
          m_trueT = std::numeric_limits<float>::quiet_NaN();
          m_trueR = std::numeric_limits<float>::quiet_NaN();

          m_residualRPhi = std::numeric_limits<float>::quiet_NaN();
          m_residualZ = std::numeric_limits<float>::quiet_NaN();
          m_residualDCA = std::numeric_limits<float>::quiet_NaN();
        } else {
          const auto& sl2 = sp.sourceLinks()[1].get<IndexSourceLink>();

          const Acts::Surface* surface2 =
              m_cfg.trackingGeometry->findSurface(sl2.geometryId());
          const auto hitRange2 =
              makeRange(measurementSimHitsMap.equal_range(m_measurementId2));
          const auto [trueLocal2, truePosition2, trueDir2] = averageSimHits(
              ctx.geoContext, *surface2, simHits, hitRange2, logger());

          const Acts::Vector4 closestPoint = closestPointOnLine(
              truePosition1, truePosition2, recoPosition.head<3>());
          const Acts::Vector4 referencePosition =
              surface1->isOnSurface(ctx.geoContext, truePosition1.head<3>(),
                                    trueDir1)
                  ? truePosition1
                  : (surface2->isOnSurface(ctx.geoContext,
                                           truePosition2.head<3>(), trueDir2)
                         ? truePosition2
                         : closestPoint);

          m_trueX = referencePosition.x() / Acts::UnitConstants::mm;
          m_trueY = referencePosition.y() / Acts::UnitConstants::mm;
          m_trueZ = referencePosition.z() / Acts::UnitConstants::mm;
          m_trueT = referencePosition.w() / Acts::UnitConstants::mm;
          m_trueR = Acts::VectorHelpers::perp(referencePosition) /
                    Acts::UnitConstants::mm;

          m_residualRPhi =
              m_trueR * Acts::detail::difference_periodic(
                            Acts::VectorHelpers::phi(recoPosition),
                            Acts::VectorHelpers::phi(referencePosition),
                            2 * std::numbers::pi);
          m_residualZ = m_trueZ - recoPosition.z() / Acts::UnitConstants::mm;
          m_residualDCA =
              (closestPoint.head<3>() - recoPosition.head<3>()).norm() /
              Acts::UnitConstants::mm;
        }
      }
    }

    // Fill the tree
    m_outputTree->Fill();
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
