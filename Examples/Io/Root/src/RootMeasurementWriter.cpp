// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"

#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <cstddef>
#include <ios>
#include <stdexcept>
#include <utility>
#include <variant>

#include <TFile.h>

namespace Acts {
class Surface;
}  // namespace Acts

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

  if (m_cfg.surfaceByIdentifier.empty()) {
    throw std::invalid_argument("Missing Surface-GeoID association map");
  }
  // Setup ROOT File
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }

  m_outputFile->cd();

  // Analyze the smearers
  std::vector<
      std::pair<Acts::GeometryIdentifier, std::unique_ptr<DigitizationTree>>>
      dTrees;
  if (!m_cfg.boundIndices.empty()) {
    ACTS_DEBUG("Bound indices are declared, preparing trees.");
    for (std::size_t ikv = 0; ikv < m_cfg.boundIndices.size(); ++ikv) {
      auto geoID = m_cfg.boundIndices.idAt(ikv);
      auto bIndices = m_cfg.boundIndices.valueAt(ikv);
      auto dTree = std::make_unique<DigitizationTree>(geoID);
      for (const auto& bIndex : bIndices) {
        ACTS_VERBOSE("- setup branch for index: " << bIndex);
        dTree->setupBoundRecBranch(bIndex);
      }
      if (!m_cfg.inputClusters.empty()) {
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
  if (!m_cfg.inputClusters.empty()) {
    clusters = m_inputClusters(ctx);
  }

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    const auto& meas = measurements[hitIdx];

    std::visit(
        [&](const auto& m) {
          Acts::GeometryIdentifier geoId =
              m.sourceLink().template get<IndexSourceLink>().geometryId();
          // find the corresponding surface
          auto surfaceItr = m_cfg.surfaceByIdentifier.find(geoId);
          if (surfaceItr == m_cfg.surfaceByIdentifier.end()) {
            return;
          }
          const Acts::Surface& surface = *(surfaceItr->second);
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
          auto [local, pos4, dir] = averageSimHits(ctx.geoContext, surface,
                                                   simHits, indices, logger());
          Acts::RotationMatrix3 rot =
              surface
                  .referenceFrame(ctx.geoContext, pos4.segment<3>(Acts::ePos0),
                                  dir)
                  .inverse();
          std::pair<double, double> angles =
              Acts::VectorHelpers::incidentAngles(dir, rot);
          dTree->fillTruthParameters(local, pos4, dir, angles);
          dTree->fillBoundMeasurement(m);
          if (!clusters.empty()) {
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
