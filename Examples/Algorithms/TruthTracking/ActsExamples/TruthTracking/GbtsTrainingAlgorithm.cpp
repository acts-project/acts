// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/GbtsTrainingAlgorithm.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsPlugins/Json/GbtsConfigJsonConverter.hpp"

#include <algorithm>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ActsExamples {

GbtsTrainingAlgorithm::GbtsTrainingAlgorithm(
    const Config& config, std::unique_ptr<const Acts::Logger> inputLogger)
    : IAlgorithm("GbtsTrainingAlgorithm", std::move(inputLogger)),
      m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }

  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurements collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing input simulated hits collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing input simulated hits measurements map");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputParticleMeasurementsMap.initialize(m_cfg.inputParticleMeasurementsMap);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);

  ACTS_INFO("LayerConnectionTool chosen");

  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  // the layers and the surfaces they are made of
  const auto layers = Acts::Experimental::readGbtsLayers(m_cfg.geometryFileDir);
  std::vector<Acts::GeometryHierarchyMap<
      Acts::Experimental::GbtsExperimentLayerId>::InputElement>
      surfaceLayers;
  for (const auto& layer : layers) {
    for (const Acts::GeometryIdentifier& surface : layer.surfaces) {
      surfaceLayers.emplace_back(surface, layer.id);
    }
  }
  m_surfaceLayers =
      Acts::GeometryHierarchyMap<Acts::Experimental::GbtsExperimentLayerId>(
          std::move(surfaceLayers));

  // the sensitive surfaces of every layer
  std::map<Acts::Experimental::GbtsExperimentLayerId,
           std::vector<const Acts::Surface*>>
      layerSurfaces;
  m_cfg.trackingGeometry->visitSurfaces([&](const Acts::Surface* surface) {
    // the entry of the module or, without one, the one of its whole layer
    const auto layer = m_surfaceLayers.find(surface->geometryId());

    if (layer == m_surfaceLayers.end()) {
      ACTS_DEBUG("No GBTS layer for volume: "
                 << surface->geometryId().volume()
                 << " Layer: " << surface->geometryId().layer()
                 << " Surface: " << surface->geometryId().sensitive());
      return;
    }

    layerSurfaces[*layer].push_back(surface);
  });

  // the symmetrization finds the mirrored layer through the r and z extent of
  // every layer, measured by a proto layer of its surfaces, which takes the
  // closest approach of a surface to the beam line and the thickness of a
  // sensitive surface into account
  const auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  auto& detectorGeometry = m_cfg.gbtsLayerConnectionToolConfig.detectorGeometry;
  detectorGeometry.clear();
  for (const auto& layer : layers) {
    const auto surfaces = layerSurfaces.find(layer.id);
    if (surfaces == layerSurfaces.end()) {
      ACTS_WARNING("No surface of GBTS layer " << layer.id
                                               << " is in the geometry");
      continue;
    }
    using enum Acts::AxisDirection;
    const Acts::ProtoLayer protoLayer(gctx, surfaces->second);
    detectorGeometry.push_back(
        {.minR = static_cast<float>(protoLayer.min(AxisR)),
         .maxR = static_cast<float>(protoLayer.max(AxisR)),
         .minZ = static_cast<float>(protoLayer.min(AxisZ)),
         .maxZ = static_cast<float>(protoLayer.max(AxisZ)),
         .gbtsId = layer.id});
    ACTS_DEBUG("GBTS layer " << layer.id << ": " << protoLayer.extent);
  }

  m_layerConnectionTool.emplace(
      m_cfg.gbtsLayerConnectionToolConfig,
      this->logger().cloneWithSuffix("GbtsLayerConnectionTool"));
}

ProcessCode GbtsTrainingAlgorithm::finalize() {
  const auto layerTable = m_layerConnectionTool->createConnectionTable();

  std::vector<Acts::Experimental::GbtsLayerConnection> connections;
  for (const auto& layerPair : layerTable) {
    // swap order as we want outward -> inward ordering
    connections.push_back({.src = layerPair.second, .dst = layerPair.first});
  }
  Acts::Experimental::writeGbtsConnections(m_cfg.outputFileDir, connections);

  return ProcessCode::SUCCESS;
}

ProcessCode GbtsTrainingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // prepare input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);
  const auto& measurementsIn = m_inputMeasurements(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& measurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  ACTS_VERBOSE("analysis hit information for " << particles.size()
                                               << " particles");

  for (const auto& [i, particle] : Acts::enumerate(particles)) {
    // find the corresponding hits for this particle
    const auto& measurements =
        makeRange(particleMeasurementsMap.equal_range(particle.particleId()));
    ACTS_VERBOSE(measurements.size()
                 << " measurements for particle " << particle);

    // the time and the GBTS layer of every hit on one of the layers
    std::vector<std::pair<double, Acts::Experimental::GbtsExperimentLayerId>>
        hits;

    hits.reserve(measurements.size());

    for (const auto& [barcode, index] : measurements) {
      ConstVariableBoundMeasurementProxy measurement =
          measurementsIn.getMeasurement(index);

      ACTS_VERBOSE("   - Measurement " << index << " with barcode " << barcode
                                       << " at " << measurement.geometryId());

      const auto simHitMapIt = measurementSimHitsMap.find(index);
      if (simHitMapIt == measurementSimHitsMap.end()) {
        ACTS_WARNING("No sim hit found for measurement index " << index);
        continue;
      }

      const auto simHitIndex = simHitMapIt->second;

      const auto simHitIt = simHits.nth(simHitIndex);
      if (simHitIt == simHits.end()) {
        ACTS_WARNING("No sim hit found for sim hit index "
                     << simHitIndex << " from measurement " << index);
        continue;
      }

      const auto layer = m_surfaceLayers.find(measurement.geometryId());
      if (layer == m_surfaceLayers.end()) {
        ACTS_DEBUG("No GBTS layer for volume: "
                   << measurement.geometryId().volume()
                   << " Layer: " << measurement.geometryId().layer()
                   << " Surface: " << measurement.geometryId().sensitive());
        continue;
      }

      hits.emplace_back(simHitIt->time(), *layer);
    }

    // the layers in the order the particle passed them
    std::ranges::sort(hits, {}, [](const auto& hit) { return hit.first; });

    std::vector<Acts::Experimental::GbtsExperimentLayerId> layers;

    layers.reserve(hits.size());

    for (const auto& [time, layer] : hits) {
      layers.push_back(layer);
    }

    {
      std::lock_guard<std::mutex> lock(m_gbtsLayerConnectionToolMutex);
      m_layerConnectionTool->addTrack(layers);
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
