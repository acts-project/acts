// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/GraphBasedSeedingAlgorithm.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/GbtsBinConnection.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsTauLookupTable.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsPlugins/Json/GbtsConfigJsonConverter.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace ActsExamples {

GraphBasedSeedingAlgorithm::GraphBasedSeedingAlgorithm(
    const Config &cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("GraphBasedSeedingAlgorithm", std::move(logger)), m_cfg(cfg) {
  // initialise the space point, seed and cluster handles
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_inputClusters.initialize(m_cfg.inputClusters);

  // parse the mapping file and turn into map
  m_actsGbtsMap = makeGbtsLayerMap(m_cfg.layerMappingFile);

  // read which eta bins may be connected
  auto connections =
      Acts::Experimental::readGbtsConnections(m_cfg.connectorInputFile);

  // keep the connections between layers of the seeded technology
  const Acts::Experimental::GbtsLayerTechnology seededTechnology =
      m_cfg.useStripConnections
          ? Acts::Experimental::GbtsLayerTechnology::Strip
          : Acts::Experimental::GbtsLayerTechnology::Pixel;
  std::map<Acts::Experimental::GbtsExperimentLayerId,
           Acts::Experimental::GbtsLayerTechnology>
      layerTechnology;
  for (const GbtsLayerInfo &layerInfo : m_actsGbtsMap) {
    layerTechnology.emplace(layerInfo.layerId, layerInfo.technology);
  }
  const auto isSeededTechnology =
      [&](Acts::Experimental::GbtsExperimentLayerId id) {
        const auto technology = layerTechnology.find(id);
        return technology != layerTechnology.end() &&
               technology->second == seededTechnology;
      };
  std::erase_if(connections,
                [&](const Acts::Experimental::GbtsBinConnection &connection) {
                  return !isSeededTechnology(connection.src.layer) ||
                         !isSeededTechnology(connection.dst.layer);
                });

  // the cluster width cuts are the only user of the tau lookup table
  if (m_cfg.seedFinderConfig.useClusterWidthCuts) {
    m_cfg.seedFinderConfig.tauLookupTable =
        Acts::Experimental::readGbtsTauLookupTable(m_cfg.lutInputFile);
  }

  // create the TrigInDetSiLayers (Logical Layers),
  // as well as a map that tracks there index in m_layerGeometry
  const auto layerGeometry = makeGbtsLayerDescriptions(
      *m_cfg.trackingGeometry, m_actsGbtsMap,
      Acts::GeometryContext::dangerouslyDefaultConstruct(), m_cfg.useLayerTilt,
      m_cfg.fillModuleCsv, this->logger());

  // initialise the object that holds all the geometry information needed for
  // the algorithm
  m_geometry = std::make_shared<Acts::Experimental::GbtsGeometry>(
      layerGeometry, connections, m_cfg.etaBinWidth, m_cfg.gbtsZ0Range,
      this->logger());

  // ROI file:Defines what region in detector we are interested in, currently
  // set to entire detector
  // for pixel seeding, roi z bounds are used
  m_internalRoi.emplace(-4.5, 4.5, m_cfg.collisionRegionMin,
                        m_cfg.collisionRegionMax);

  // the RoI owns the luminous region, so it overrides what came in
  m_cfg.graphConfig.maxZ0 = m_internalRoi->zMax();
  m_cfg.graphConfig.minZ0 = m_internalRoi->zMin();

  m_gbtsGraphBuilder.emplace(
      m_cfg.graphConfig, m_geometry,
      this->logger().cloneWithSuffix("GbtsGraphBuilder"));

  m_finder.emplace(Acts::Experimental::GraphBasedTrackSeeder::DerivedConfig(
                       m_cfg.seedFinderConfig),
                   m_geometry, this->logger().cloneWithSuffix("GbtsFinder"));

  m_filter = Acts::Experimental::GbtsTrackingFilter(
      m_cfg.trackingFilterConfig, m_geometry,
      this->logger().cloneWithSuffix("GbtsFilter"));

  printConfig();
}

ProcessCode GraphBasedSeedingAlgorithm::execute(
    const AlgorithmContext &ctx) const {
  // initialise input space points from handle and define new container
  const SpacePointContainer &spacePoints = m_inputSpacePoints(ctx);

  const Acts::Experimental::GraphBasedTrackSeeder::Options options{
      .bFieldInZ = m_cfg.bFieldInZ};

  // The node storage is filled straight from the input space points. It takes
  // plain scalars, so no intermediate space point container is needed and the
  // seeds come back indexed into the input container directly.
  Acts::Experimental::GbtsNodeStorage nodeStorage = m_finder->makeNodeStorage();

  std::uint32_t nUnmapped = 0;

  for (const auto &spacePoint : spacePoints) {
    const std::optional<Acts::Experimental::GbtsLayerIndex> layerIndex =
        gbtsLayerIndex(spacePoint);
    if (!layerIndex.has_value()) {
      ++nUnmapped;
      continue;
    }

    // Cluster width and local position are not available in the examples
    // framework, so the machine learning features stay switched off.
    nodeStorage.insert(spacePoint.index(), spacePoint.x(), spacePoint.y(),
                       spacePoint.z(), *layerIndex);
  }

  nodeStorage.finalize();

  ACTS_VERBOSE("Loaded " << nodeStorage.numberOfNodes() << " graph nodes, "
                         << nUnmapped << " space points not in the GBTS map");

  Acts::SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);

  // create the seeds

  m_finder->createSeeds(nodeStorage, m_internalRoi.value(), *m_gbtsGraphBuilder,
                        *m_filter, options, seeds);

  m_outputSeeds(ctx, std::move(seeds));

  return ProcessCode::SUCCESS;
}

std::optional<Acts::Experimental::GbtsLayerIndex>
GraphBasedSeedingAlgorithm::gbtsLayerIndex(
    const ConstSpacePointProxy &spacePoint) const {
  const auto &sourceLink = spacePoint.sourceLinks();

  if (sourceLink.empty()) {
    ACTS_WARNING("warning source link vector is empty");
    return std::nullopt;
  }

  const auto &indexSourceLink = sourceLink.front().get<IndexSourceLink>();

  const Acts::GeometryIdentifier geoId = indexSourceLink.geometryId();

  // the entry of the module or, without one, the one of its whole layer
  const auto find = m_actsGbtsMap.find(geoId);

  // a space point off the GBTS layers takes no part in the seeding
  if (find == m_actsGbtsMap.end()) {
    ACTS_DEBUG("No GBTS layer for volume: "
               << geoId.volume() << " Layer: " << geoId.layer()
               << " Surface: " << geoId.sensitive());
    return std::nullopt;
  }

  return m_geometry->layerIndex(find->layerId);
}

void GraphBasedSeedingAlgorithm::printConfig() const {
  ACTS_DEBUG("===== GraphBasedSeedingAlgorithm =====");
  ACTS_DEBUG("layerMappingFile: " << m_cfg.layerMappingFile);
  ACTS_DEBUG("connectorInputFile: " << m_cfg.connectorInputFile);
  ACTS_DEBUG("lutInputFile: " << m_cfg.lutInputFile);
  ACTS_DEBUG("etaBinWidth: " << m_cfg.etaBinWidth);
  ACTS_DEBUG("collisionRegionMin: " << m_cfg.collisionRegionMin);
  ACTS_DEBUG("collisionRegionMax: " << m_cfg.collisionRegionMax);
  ACTS_DEBUG("useStripConnections: " << m_cfg.useStripConnections);
  ACTS_DEBUG("===== GraphBasedTrackSeeder =====");
  const auto &cfg1 = m_cfg.seedFinderConfig;
  ACTS_DEBUG("nMaxPhiSlice: " << cfg1.nMaxPhiSlice);
  ACTS_DEBUG("useClusterWidthCuts: " << cfg1.useClusterWidthCuts);
  ACTS_DEBUG("edgeMaskMinEta: " << cfg1.edgeMaskMinEta);
  ACTS_DEBUG("hitShareThreshold: " << cfg1.hitShareThreshold);
  ACTS_DEBUG("maxEndcapClusterWidth: " << cfg1.maxEndcapClusterWidth);
  ACTS_DEBUG("maxSeedSplitEta: " << cfg1.maxSeedSplitEta);
  ACTS_DEBUG("maxInvRadDiff: " << cfg1.maxInvRadDiff);
  ACTS_DEBUG("ccaMaxIterations: " << cfg1.ccaMaxIterations);
  ACTS_DEBUG("minSeedLevel: " << static_cast<std::uint32_t>(cfg1.minSeedLevel));
  ACTS_DEBUG("addTriplets: " << cfg1.addTriplets);
  ACTS_DEBUG("maxAbsEtaAddTriplets: " << cfg1.maxAbsEtaAddTriplets);
  ACTS_DEBUG("=====GbtsGraphBuilder=====");
  const auto &cfg2 = m_cfg.graphConfig;
  ACTS_DEBUG("matchBeforeCreate: " << cfg2.matchBeforeCreate);
  ACTS_DEBUG("tauRatioCut: " << cfg2.tauRatioCut);
  ACTS_DEBUG("tauRatioPrecut: " << cfg2.tauRatioPrecut);
  ACTS_DEBUG("minPt: " << cfg2.minPt);
  ACTS_DEBUG("useEtaBinning: " << cfg2.useEtaBinning);
  ACTS_DEBUG("doubletFilterRZ: " << cfg2.doubletFilterRZ);
  ACTS_DEBUG("maxEdgesPerSP: " << cfg2.maxEdgesPerSP);
  ACTS_DEBUG("minDeltaRadius: " << cfg2.minDeltaRadius);
  ACTS_DEBUG("validateTriplets: " << cfg2.validateTriplets);
  ACTS_DEBUG("useAdaptiveCuts: " << cfg2.useAdaptiveCuts);
  ACTS_DEBUG("tauRatioCorr: " << cfg2.tauRatioCorr);
  ACTS_DEBUG("d0Max: " << cfg2.d0Max);
  ACTS_DEBUG("cutDPhiMax: " << cfg2.cutDPhiMax);
  ACTS_DEBUG("cutDCurvMax: " << cfg2.cutDCurvMax);
  ACTS_DEBUG("minZ0: " << cfg2.minZ0);
  ACTS_DEBUG("maxZ0: " << cfg2.maxZ0);
  ACTS_DEBUG("maxOuterRadius: " << cfg2.maxOuterRadius);
  ACTS_DEBUG("===== GbtsTrackFilter =====");
  const auto &cfg3 = m_cfg.trackingFilterConfig;
  ACTS_DEBUG("sigmaMS: " << cfg3.sigmaMS);
  ACTS_DEBUG("radLen: " << cfg3.radLen);
  ACTS_DEBUG("sigmaX: " << cfg3.sigmaX);
  ACTS_DEBUG("sigmaY: " << cfg3.sigmaY);
  ACTS_DEBUG("weightX: " << cfg3.weightX);
  ACTS_DEBUG("weightY: " << cfg3.weightY);
  ACTS_DEBUG("maxDChi2X: " << cfg3.maxDChi2X);
  ACTS_DEBUG("maxDChi2Y: " << cfg3.maxDChi2Y);
  ACTS_DEBUG("addHit: " << cfg3.addHit);
  ACTS_DEBUG("maxCurvature: " << cfg3.maxCurvature);
  ACTS_DEBUG("maxZ0: " << cfg3.maxZ0);
  ACTS_DEBUG("================================");
}

}  // namespace ActsExamples
