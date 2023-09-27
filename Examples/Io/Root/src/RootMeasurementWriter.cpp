// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <ios>
#include <optional>
#include <stdexcept>

#include <TFile.h>
#include <TString.h>

ActsExamples::RootMeasurementWriter::RootMeasurementWriter(
    const ActsExamples::RootMeasurementWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputMeasurements, "RootMeasurementWriter", level),
      m_cfg(config) {
  // Input container for measurements is already checked by base constructor
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map input collection");
  }

  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);

  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  // Setup ROOT File
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
  }

  m_outputFile->cd();

  // Analyze the smearers
  std::vector<
      std::pair<Acts::GeometryIdentifier, std::unique_ptr<DigitizationTree>>>
      dTrees;
  if (not m_cfg.boundIndices.empty()) {
    ACTS_DEBUG("Bound indices are declared, preparing trees.");
    for (size_t ikv = 0; ikv < m_cfg.boundIndices.size(); ++ikv) {
      auto geoID = m_cfg.boundIndices.idAt(ikv);
      auto bIndices = m_cfg.boundIndices.valueAt(ikv);
      auto dTree = std::make_unique<DigitizationTree>(geoID);
      for (const auto& bIndex : bIndices) {
        ACTS_VERBOSE("- setup branch for index: " << bIndex);
        dTree->setupBoundRecBranch(bIndex);
      }
      if (not m_cfg.inputClusters.empty()) {
        dTree->setupClusterBranch(bIndices);
      }
      dTrees.push_back({geoID, std::move(dTree)});
    }
  } else {
    ACTS_DEBUG("Bound indices are not declared, no reco setup.")
  }

  m_outputTrees = Acts::GeometryHierarchyMap<std::unique_ptr<DigitizationTree>>(
      std::move(dTrees));
}

ActsExamples::RootMeasurementWriter::~RootMeasurementWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootMeasurementWriter::finalize() {
  /// Close the file if it's yours
  m_outputFile->cd();
  for (auto dTree = m_outputTrees.begin(); dTree != m_outputTrees.end();
       ++dTree) {
    (*dTree)->tree->Write();
  }
  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootMeasurementWriter::writeT(
    const AlgorithmContext& ctx, const MeasurementContainer& measurements) {
  const auto& simHits = m_inputSimHits(ctx);
  const auto& hitSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  ClusterContainer clusters;
  if (not m_cfg.inputClusters.empty()) {
    clusters = m_inputClusters(ctx);
  }

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    const auto& meas = measurements[hitIdx];

    std::visit(
        [&](const auto& m) {
          Acts::GeometryIdentifier geoId = m.sourceLink().geometryId();
          // find the corresponding surface
          const Acts::Surface* surfacePtr =
              m_cfg.trackingGeometry->findSurface(geoId);
          if (not surfacePtr) {
            return;
          }
          const Acts::Surface& surface = *surfacePtr;
          // find the corresponding output tree
          auto dTreeItr = m_outputTrees.find(geoId);
          if (dTreeItr == m_outputTrees.end()) {
            return;
          }
          auto& dTree = *dTreeItr;

          // Fill the identification
          dTree->fillIdentification(ctx.eventNumber, geoId);

          // Find the contributing simulated hits
          auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
          // Use average truth in the case of multiple contributing sim hits
          auto [local, pos4, dir] =
              averageSimHits(ctx.geoContext, surface, simHits, indices);
          dTree->fillTruthParameters(local, pos4, dir);
          dTree->fillBoundMeasurement(m);
          if (not clusters.empty()) {
            const auto& c = clusters[hitIdx];
            dTree->fillCluster(c);
          }
          dTree->tree->Fill();
          if (dTree->chValue != nullptr) {
            dTree->chValue->clear();
          }
          if (dTree->chId[0] != nullptr) {
            dTree->chId[0]->clear();
          }
          if (dTree->chId[1] != nullptr) {
            dTree->chId[1]->clear();
          }
        },
        meas);
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
