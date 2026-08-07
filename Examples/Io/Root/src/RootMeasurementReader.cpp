// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMeasurementReader.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/Digitization/Segmentizer.hpp"
#include "ActsPlugins/Root/RootMeasurementIo.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

#include <TChain.h>

namespace ActsExamples {

RootMeasurementReader::RootMeasurementReader(
    const RootMeasurementReader::Config& config, Acts::Logging::Level level)
    : m_cfg(config), m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement output collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing input filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputMeasurementSubset.maybeInitialize(m_cfg.outputMeasurementSubset);
  m_outputMeasurementParticlesMap.maybeInitialize(
      m_cfg.outputMeasurementParticlesMap);
  m_outputClusters.maybeInitialize(m_cfg.outputClusters);

  m_inputChain = std::make_unique<TChain>(m_cfg.treeName.c_str());
  m_inputChain->Add(m_cfg.filePath.c_str());
  m_inputChain->LoadTree(0);
  ACTS_DEBUG("Adding File " << m_cfg.filePath << " to tree '" << m_cfg.treeName
                            << "'.");

  ActsPlugins::RootMeasurementIo::Config ioCfg{m_cfg.boundIndices,
                                               m_cfg.clusterIndices};
  m_measurementIo = std::make_unique<ActsPlugins::RootMeasurementIo>(ioCfg);
  m_measurementIo->connectForRead(*m_inputChain);

  m_inputChain->SetBranchAddress("particles_vertex_primary",
                                 &m_particleVertexPrimary.get());
  m_inputChain->SetBranchAddress("particles_vertex_secondary",
                                 &m_particleVertexSecondary.get());
  m_inputChain->SetBranchAddress("particles_particle",
                                 &m_particleParticle.get());
  m_inputChain->SetBranchAddress("particles_generation",
                                 &m_particleGeneration.get());
  m_inputChain->SetBranchAddress("particles_sub_particle",
                                 &m_particleSubParticle.get());

  // Because each measurement is stored in a single entry in the root file, we
  // need to scan the file first for the positions of the events in the file
  // in order to efficiently read the events later on.

  // Disable all branches and only enable event_nr for a first scan of the
  // file
  m_inputChain->SetBranchStatus("*", false);
  m_inputChain->SetBranchStatus("event_nr", true);

  auto nEntries = static_cast<std::size_t>(m_inputChain->GetEntriesFast());
  if (nEntries == 0) {
    throw std::runtime_error("Did not find any entries in input file");
  }

  m_inputChain->GetEntry(0);
  m_eventMap.push_back(
      {static_cast<std::uint32_t>(m_measurementIo->eventNr()), 0ul, 0ul});

  for (auto i = 1ul; i < nEntries; ++i) {
    m_inputChain->GetEntry(i);
    const auto evtNr = static_cast<std::uint32_t>(m_measurementIo->eventNr());

    if (evtNr != std::get<0>(m_eventMap.back())) {
      std::get<2>(m_eventMap.back()) = i;
      m_eventMap.push_back({evtNr, i, i});
    }
  }
  std::get<2>(m_eventMap.back()) = nEntries;

  std::ranges::sort(m_eventMap, {},
                    [](const auto& m) { return std::get<0>(m); });

  // Re-enable all branches
  m_inputChain->SetBranchStatus("*", true);
  ACTS_DEBUG("Event range: " << availableEvents().first << " - "
                             << availableEvents().second);
}

RootMeasurementReader::~RootMeasurementReader() = default;

std::pair<std::size_t, std::size_t> RootMeasurementReader::availableEvents()
    const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode RootMeasurementReader::read(const AlgorithmContext& ctx) {
  auto it = std::ranges::find_if(m_eventMap, [&](const auto& a) {
    return std::get<0>(a) == ctx.eventNumber;
  });

  MeasurementContainer measurements;
  MeasurementParticlesMap measurementParticlesMap;
  ClusterContainer clusters;

  if (it == m_eventMap.end()) {
    ACTS_DEBUG("Reading empty event: " << ctx.eventNumber);
  } else {
    // lock the mutex
    std::lock_guard<std::mutex> lock(m_readMutex);

    ACTS_DEBUG("Reading event: " << std::get<0>(*it)
                                 << " stored in entries: " << std::get<1>(*it)
                                 << " - " << std::get<2>(*it));

    for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); ++entry) {
      m_inputChain->GetEntry(entry);

      const Acts::GeometryIdentifier geoId = m_measurementIo->geometryId();
      auto [values, variances, subspaceIndex] =
          m_measurementIo->boundMeasurement();

      DigitizedParameters dParams;
      dParams.values = std::move(values);
      dParams.variances = std::move(variances);
      for (auto idx : subspaceIndex) {
        dParams.indices.push_back(static_cast<Acts::BoundIndices>(idx));
      }

      auto measurement = createMeasurement(measurements, geoId, dParams);
      const auto measurementIndex = static_cast<Index>(measurement.index());

      if (m_outputMeasurementParticlesMap.isInitialized()) {
        for (std::size_t i = 0; i < m_particleVertexPrimary->size(); ++i) {
          auto barcode =
              SimBarcode()
                  .withVertexPrimary(static_cast<SimBarcode::PrimaryVertexId>(
                      m_particleVertexPrimary->at(i)))
                  .withVertexSecondary(
                      static_cast<SimBarcode::SecondaryVertexId>(
                          m_particleVertexSecondary->at(i)))
                  .withParticle(static_cast<SimBarcode::ParticleId>(
                      m_particleParticle->at(i)))
                  .withGeneration(static_cast<SimBarcode::GenerationId>(
                      m_particleGeneration->at(i)))
                  .withSubParticle(static_cast<SimBarcode::SubParticleId>(
                      m_particleSubParticle->at(i)));
          measurementParticlesMap.emplace(measurementIndex, barcode);
        }
      }

      if (m_outputClusters.isInitialized()) {
        Cluster cluster;
        auto size = m_measurementIo->clusterSize();
        cluster.sizeLoc0 = static_cast<std::size_t>(size[0]);
        cluster.sizeLoc1 = static_cast<std::size_t>(size[1]);
        cluster.globalPosition = m_measurementIo->globalPosition();
        for (auto [ch0, ch1, chValue] : m_measurementIo->clusterChannels()) {
          ActsFatras::Segmentizer::Segment2D dummySegment = {
              Acts::Vector2::Zero(), Acts::Vector2::Zero()};
          ActsFatras::Segmentizer::Bin2D bin{static_cast<unsigned int>(ch0),
                                             static_cast<unsigned int>(ch1)};
          cluster.channels.emplace_back(bin, dummySegment, chValue);
        }
        clusters.push_back(std::move(cluster));
      }
    }
  }

  ACTS_DEBUG("Read " << measurements.size() << " measurements for event "
                     << ctx.eventNumber);

  const auto& storedMeasurements =
      m_outputMeasurements(ctx, std::move(measurements));

  if (m_outputMeasurementSubset.isInitialized()) {
    // Build the full subset: all measurements, indices in original space.
    std::vector<MeasurementContainer::Index> allIndices(
        storedMeasurements.size());
    std::iota(allIndices.begin(), allIndices.end(), Index{0});
    m_outputMeasurementSubset(
        ctx, MeasurementSubset(storedMeasurements, std::move(allIndices)));
  }
  if (m_outputMeasurementParticlesMap.isInitialized()) {
    m_outputMeasurementParticlesMap(ctx, std::move(measurementParticlesMap));
  }
  if (m_outputClusters.isInitialized()) {
    m_outputClusters(ctx, std::move(clusters));
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
