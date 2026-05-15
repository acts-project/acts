// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"

#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsPlugins/Root/RootMeasurementIo.hpp"

#include <ios>
#include <memory>
#include <stdexcept>
#include <utility>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {

namespace {

std::tuple<std::vector<double>, std::vector<double>, std::vector<unsigned int>>
prepareBoundMeasurement(const ConstVariableBoundMeasurementProxy& m) {
  std::vector<double> measurements = {};
  std::vector<double> variances = {};
  std::vector<unsigned int> subspaceIndex = {};

  for (unsigned int i = 0; i < m.size(); ++i) {
    auto ib = m.subspaceIndexVector()[i];

    measurements.push_back(m.parameters()[i]);
    variances.push_back(m.covariance()(i, i));
    subspaceIndex.push_back(static_cast<unsigned int>(ib));
  }

  return {measurements, variances, subspaceIndex};
}

std::vector<std::tuple<int, int, float>> prepareChannels(const Cluster& c) {
  std::vector<std::tuple<int, int, float>> channels = {};
  for (auto ch : c.channels) {
    channels.emplace_back(static_cast<int>(ch.bin[0]),
                          static_cast<int>(ch.bin[1]),
                          static_cast<float>(ch.activation));
  }
  return channels;
}

}  // namespace

RootMeasurementWriter::RootMeasurementWriter(
    const RootMeasurementWriter::Config& config, Acts::Logging::Level level)
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

  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);

  if (m_cfg.surfaceByIdentifier.empty()) {
    throw std::invalid_argument("Missing Surface-GeoID association map");
  }
  // Setup ROOT File
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }

  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), "Measurements");
  m_outputTree->Branch("particles_vertex_primary", &m_particleVertexPrimary);
  m_outputTree->Branch("particles_vertex_secondary",
                       &m_particleVertexSecondary);
  m_outputTree->Branch("particles_particle", &m_particleParticle);
  m_outputTree->Branch("particles_generation", &m_particleGeneration);
  m_outputTree->Branch("particles_sub_particle", &m_particleSubParticle);

  ActsPlugins::RootMeasurementIo::Config treeCfg{m_cfg.boundIndices,
                                                 m_cfg.clusterIndices};
  m_measurementIo = std::make_unique<ActsPlugins::RootMeasurementIo>(treeCfg);
  m_measurementIo->connectForWrite(*m_outputTree);
}

RootMeasurementWriter::~RootMeasurementWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootMeasurementWriter::finalize() {
  /// Close the file if it's yours
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

ProcessCode RootMeasurementWriter::writeT(
    const AlgorithmContext& ctx, const MeasurementContainer& measurements) {
  const auto& simHits = m_inputSimHits(ctx);
  const auto& hitSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  const ClusterContainer* clusters = nullptr;
  if (!m_cfg.inputClusters.empty()) {
    clusters = &m_inputClusters(ctx);
  }

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    const ConstVariableBoundMeasurementProxy meas =
        measurements.getMeasurement(hitIdx);

    Acts::GeometryIdentifier geoId = meas.geometryId();
    // find the corresponding surface
    auto surfaceItr = m_cfg.surfaceByIdentifier.find(geoId);
    if (surfaceItr == m_cfg.surfaceByIdentifier.end()) {
      continue;
    }
    const Acts::Surface& surface = *(surfaceItr->second);

    // Fill the identification
    m_measurementIo->fillIdentification(static_cast<int>(ctx.eventNumber),
                                        geoId);

    // Find the contributing simulated hits
    auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
    // Use average truth in the case of multiple contributing sim hits
    auto [local, pos4, dir] =
        averageSimHits(ctx.geoContext, surface, simHits, indices, logger());
    Acts::RotationMatrix3 rot =
        surface
            .referenceFrame(ctx.geoContext, pos4.segment<3>(Acts::ePos0), dir)
            .inverse();
    std::pair<double, double> angles =
        Acts::VectorHelpers::incidentAngles(dir, rot);
    for (auto [_, i] : indices) {
      const auto barcode = simHits.nth(i)->particleId();
      m_particleVertexPrimary.push_back(barcode.vertexPrimary());
      m_particleVertexSecondary.push_back(barcode.vertexSecondary());
      m_particleParticle.push_back(barcode.particle());
      m_particleGeneration.push_back(barcode.generation());
      m_particleSubParticle.push_back(barcode.subParticle());
    }
    m_measurementIo->fillTruthParameters(local, pos4, dir, angles);

    // Fill the measurement parameters & clusters still
    auto [msm, vcs, ssi] = prepareBoundMeasurement(meas);
    m_measurementIo->fillBoundMeasurement(msm, vcs, ssi);

    // Fill the cluster information if available
    if (clusters != nullptr) {
      const auto& cluster = (*clusters)[hitIdx];
      m_measurementIo->fillGlobalPosition(cluster.globalPosition);
      m_measurementIo->fillCluster(prepareChannels(cluster));
    }

    m_outputTree->Fill();
    m_particleVertexPrimary.clear();
    m_particleVertexSecondary.clear();
    m_particleParticle.clear();
    m_particleGeneration.clear();
    m_particleSubParticle.clear();
    m_measurementIo->clear();
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
