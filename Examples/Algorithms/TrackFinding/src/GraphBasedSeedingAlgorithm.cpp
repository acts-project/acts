// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/GraphBasedSeedingAlgorithm.hpp"

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsPlugins/Json/GbtsConfigJsonConverter.hpp"
#include "ActsPlugins/Json/detail/JsonIo.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
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
  m_actsGbtsMap = makeActsGbtsMap();

  // read which layers may be connected
  auto connectorTable = Acts::detail::readJsonFile(m_cfg.connectorInputFile)
                            .get<Acts::Experimental::GbtsConnectionsConfig>();

  // keep the connections between layers of the seeded technology
  const Acts::Experimental::GbtsLayerTechnology seededTechnology =
      m_cfg.seedFinderConfig.useStripConnections
          ? Acts::Experimental::GbtsLayerTechnology::Strip
          : Acts::Experimental::GbtsLayerTechnology::Pixel;
  std::map<Acts::Experimental::GbtsExperimentLayerId,
           Acts::Experimental::GbtsLayerTechnology>
      layerTechnology;
  for (const auto &[actsId, gbtsId] : m_actsGbtsMap) {
    layerTechnology.emplace(gbtsId.layerId, gbtsId.technology);
  }
  const auto isSeededTechnology =
      [&](Acts::Experimental::GbtsExperimentLayerId id) {
        const auto technology = layerTechnology.find(id);
        return technology != layerTechnology.end() &&
               technology->second == seededTechnology;
      };
  std::erase_if(connectorTable.connections,
                [&](const Acts::Experimental::GbtsLayerConnection &connection) {
                  return !isSeededTechnology(connection.src) ||
                         !isSeededTechnology(connection.dst);
                });

  // the cluster width cuts are the only user of the tau lookup table
  if (m_cfg.seedFinderConfig.useClusterWidthCuts) {
    m_cfg.seedFinderConfig.tauLookupTable =
        Acts::detail::readJsonFile(m_cfg.lutInputFile)
            .at("tauLookupTable")
            .get<Acts::Experimental::detail::GbtsTauLookupTable>();
  }

  // create the TrigInDetSiLayers (Logical Layers),
  // as well as a map that tracks there index in m_layerGeometry
  const auto layerGeometry =
      layerNumbering(Acts::GeometryContext::dangerouslyDefaultConstruct());

  // option that allows for adding custom eta binning (default is at 0.2)
  const float etaBinWidth = m_cfg.etaBinWidthOverride != 0.0f
                                ? m_cfg.etaBinWidthOverride
                                : connectorTable.etaBinWidth;

  // initialise the object that holds all the geometry information needed for
  // the algorithm
  auto geometry = std::make_shared<Acts::Experimental::GbtsGeometry>(
      layerGeometry, connectorTable.connections, etaBinWidth, m_cfg.gbtsZ0Range,
      this->logger());

  resolveLayerIndices(*geometry);

  // ROI file:Defines what region in detector we are interested in, currently
  // set to entire detector
  // for pixel seeding, roi z bounds are used

  m_internalRoi.emplace(-4.5, 4.5, -150., 150.);
  m_cfg.seedFinderConfig.maxZ0 = m_internalRoi->zMax();
  m_cfg.seedFinderConfig.minZ0 = m_internalRoi->zMin();

  m_finder = Acts::Experimental::GraphBasedTrackSeeder(
      Acts::Experimental::GraphBasedTrackSeeder::DerivedConfig(
          m_cfg.seedFinderConfig),
      geometry, this->logger().cloneWithSuffix("GbtsFinder"));

  m_filter = Acts::Experimental::GbtsTrackingFilter(
      m_cfg.trackingFilterConfig, geometry,
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

  m_finder->createSeeds(nodeStorage, m_internalRoi.value(), *m_filter, options,
                        seeds);

  m_outputSeeds(ctx, std::move(seeds));

  return ProcessCode::SUCCESS;
}

std::map<GraphBasedSeedingAlgorithm::ActsIDs,
         GraphBasedSeedingAlgorithm::GbtsIDs>
GraphBasedSeedingAlgorithm::makeActsGbtsMap() const {
  std::map<ActsIDs, GbtsIDs> actsToGbtsMap;

  // one entry per surface of a layer, sensitive 0 for a whole geometry layer
  for (const Acts::Experimental::GbtsLayerConfig &layer :
       Acts::detail::readJsonFile(m_cfg.layerMappingFile)
           .at("layers")
           .get<std::vector<Acts::Experimental::GbtsLayerConfig>>()) {
    for (const Acts::GeometryIdentifier &surface : layer.surfaces) {
      const ActsIDs actsId{surface.volume() * 100 + surface.layer(),
                           surface.sensitive()};
      const GbtsIDs gbtsId{.layerId = layer.id,
                           .type = layer.type,
                           .technology = layer.technology};
      actsToGbtsMap.insert({actsId, gbtsId});
    }
  }

  return actsToGbtsMap;
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

  const auto actsVolId =
      static_cast<std::uint32_t>(indexSourceLink.geometryId().volume());
  const auto actsLayId =
      static_cast<std::uint32_t>(indexSourceLink.geometryId().layer());
  const auto actsModId =
      static_cast<std::uint32_t>(indexSourceLink.geometryId().sensitive());

  // Search for vol, lay and module=0, if doesn't esist (end) then search
  // for full thing vol*100+lay as first number in pair then 0 or mod id
  const std::uint64_t actsJointId = std::uint64_t{actsVolId} * 100 + actsLayId;

  // here the key needs to be pair of(vol*100+lay, 0)
  ActsIDs key{actsJointId, 0};
  auto find = m_actsGbtsMap.find(key);

  // if end then make new key of (vol*100+lay, modid)
  if (find == m_actsGbtsMap.end()) {
    key = ActsIDs{actsJointId, actsModId};  // mod ID
    find = m_actsGbtsMap.find(key);
  }

  // a space point off the GBTS layers takes no part in the seeding
  if (find == m_actsGbtsMap.end()) {
    ACTS_DEBUG("Key not found in Gbts map for volume id: "
               << actsVolId << " and layer id: " << actsLayId);
    return std::nullopt;
  }

  return find->second.layerIndex;
}

void GraphBasedSeedingAlgorithm::resolveLayerIndices(
    const Acts::Experimental::GbtsGeometry &geometry) {
  for (auto &[actsId, gbtsId] : m_actsGbtsMap) {
    const std::optional<Acts::Experimental::GbtsLayerIndex> index =
        geometry.layerIndex(gbtsId.layerId);

    if (!index.has_value()) {
      ACTS_WARNING("No GBTS layer for ID: " << gbtsId.layerId);
    }

    gbtsId.layerIndex = index;
  }
}

std::vector<Acts::Experimental::GbtsLayerDescription>
GraphBasedSeedingAlgorithm::layerNumbering(
    const Acts::GeometryContext &gctx) const {
  std::vector<Acts::Experimental::GbtsLayerDescription> inputVector;
  std::vector<std::size_t> countVector;

  m_cfg.trackingGeometry->visitSurfaces(
      [this, &inputVector, &countVector, &gctx](const Acts::Surface *surface) {
        Acts::GeometryIdentifier geoId = surface->geometryId();
        auto actsVolId = geoId.volume();
        auto actsLayId = geoId.layer();
        auto mod_id = geoId.sensitive();
        auto bounds_vect = surface->bounds().values();
        auto center = surface->center(gctx);

        // make bounds global
        Acts::Vector3 globalFakeMom(1, 1, 1);
        Acts::Vector2 min_bound_local =
            Acts::Vector2(bounds_vect[0], bounds_vect[1]);
        Acts::Vector2 max_bound_local =
            Acts::Vector2(bounds_vect[2], bounds_vect[3]);
        Acts::Vector3 min_bound_global =
            surface->localToGlobal(gctx, min_bound_local, globalFakeMom);
        Acts::Vector3 max_bound_global =
            surface->localToGlobal(gctx, max_bound_local, globalFakeMom);

        // checking that not wrong way round
        if (min_bound_global(0) > max_bound_global(0)) {
          min_bound_global.swap(max_bound_global);
        }

        float rc = 0.0;
        float minBound = 100000.0;
        float maxBound = -100000.0;

        // convert to Gbts ID
        auto actsJointId = actsVolId * 100 + actsLayId;
        // here the key needs to be pair of(vol*100+lay, 0)
        auto key = ActsIDs{actsJointId, 0};
        auto find = m_actsGbtsMap.find(key);

        // check to see if key exists
        if (find == m_actsGbtsMap.end()) {
          key = ActsIDs{actsJointId, mod_id};
          find = m_actsGbtsMap.find(key);
        }

        // a surface off the GBTS layers takes no part in the seeding
        if (find == m_actsGbtsMap.end()) {
          ACTS_DEBUG("Key not found in Gbts map for volume id: "
                     << actsVolId << ", layer id: " << actsLayId
                     << ", sensitive id: " << mod_id);
          return;  // skip this surface in the visitor
        }

        const Acts::Experimental::GbtsExperimentLayerId gbtsId =
            find->second.layerId;

        // a variable that says if barrrel, 0 = barrel
        Acts::Experimental::GbtsLayerType barrelEc = find->second.type;

        if (barrelEc == Acts::Experimental::GbtsLayerType::Barrel) {
          rc = std::sqrt(center(0) * center(0) +
                         center(1) * center(1));  // barrel center in r
          // bounds of z
          if (min_bound_global(2) < minBound) {
            minBound = min_bound_global(2);
          }
          if (max_bound_global(2) > maxBound) {
            maxBound = max_bound_global(2);
          }
        } else if (barrelEc == Acts::Experimental::GbtsLayerType::Endcap) {
          rc = center(2);  // not barrel center in Z
          // bounds of r
          float min = std::sqrt(min_bound_global(0) * min_bound_global(0) +
                                min_bound_global(1) * min_bound_global(1));
          float max = std::sqrt(max_bound_global(0) * max_bound_global(0) +
                                max_bound_global(1) * max_bound_global(1));
          if (min < minBound) {
            minBound = min;
          }
          if (max > maxBound) {
            maxBound = max;
          }
        } else {
          throw std::runtime_error(
              "Invalid barrel/endcap assignment for GbtsLayer");
        }

        const auto currentIndex =
            find_if(inputVector.begin(), inputVector.end(),
                    [gbtsId](auto n) { return n.id == gbtsId; });
        if (currentIndex != inputVector.end()) {  // not end so does exist
          const auto index = static_cast<std::size_t>(
              std::distance(inputVector.begin(), currentIndex));
          inputVector[index].refCoord += rc;
          inputVector[index].minBound =
              std::min(inputVector[index].minBound, minBound);
          inputVector[index].maxBound =
              std::max(inputVector[index].maxBound, maxBound);
          countVector[index] += 1;  // increase count at the index

        } else {  // end so doesn't exists
          // make new if one with Gbts ID doesn't exist:
          inputVector.push_back(Acts::Experimental::GbtsLayerDescription{
              .id = gbtsId,
              .type = barrelEc,
              .technology = find->second.technology,
              .refCoord = rc,
              .minBound = minBound,
              .maxBound = maxBound});
          // so the element exists and not divinding by 0
          countVector.push_back(1);
        }

        // add to file each time,
        // print to csv for each module, no repeats so dont need to make
        // map for averaging
        if (m_cfg.fillModuleCsv) {
          std::fstream fout;
          fout.open("ACTS_modules.csv", std::ios::out | std::ios::app);
          fout << actsVolId << ", "  // vol
               << actsLayId << ", "  // lay
               << mod_id << ", "     // module
               << gbtsId << ","      // Gbts id
               << center(2) << ", "  // z
               << std::sqrt(center(0) * center(0) + center(1) * center(1))  // r
               << "\n";
        }
      });

  for (std::size_t i = 0; i < inputVector.size(); i++) {
    inputVector[i].refCoord = inputVector[i].refCoord / countVector[i];
  }

  return inputVector;
}

void GraphBasedSeedingAlgorithm::printConfig() const {
  ACTS_DEBUG("===== GraphBasedSeedingAlgorithm =====");
  ACTS_DEBUG("layerMappingFile: " << m_cfg.layerMappingFile);
  ACTS_DEBUG("connectorInputFile: " << m_cfg.connectorInputFile);
  ACTS_DEBUG("lutInputFile: " << m_cfg.lutInputFile);
  ACTS_DEBUG("etaBinWidthOverride: " << m_cfg.etaBinWidthOverride);
  ACTS_DEBUG("===== GraphBasedTrackSeeder =====");
  const auto &cfg1 = m_cfg.seedFinderConfig;
  ACTS_DEBUG("BeamSpotCorrection: " << cfg1.beamSpotCorrection);
  ACTS_DEBUG("useStripConnections: " << cfg1.useStripConnections);
  ACTS_DEBUG("useClusterWidthCuts: " << cfg1.useClusterWidthCuts);
  ACTS_DEBUG("matchBeforeCreate: " << cfg1.matchBeforeCreate);
  ACTS_DEBUG("tauRatioCut: " << cfg1.tauRatioCut);
  ACTS_DEBUG("tauRatioPrecut: " << cfg1.tauRatioPrecut);
  ACTS_DEBUG("nMaxPhiSlice: " << cfg1.nMaxPhiSlice);
  ACTS_DEBUG("minPt: " << cfg1.minPt);
  ACTS_DEBUG("useEtaBinning: " << cfg1.useEtaBinning);
  ACTS_DEBUG("doubletFilterRZ: " << cfg1.doubletFilterRZ);
  ACTS_DEBUG("nMaxEdges: " << cfg1.nMaxEdges);
  ACTS_DEBUG("minDeltaRadius: " << cfg1.minDeltaRadius);
  ACTS_DEBUG("edgeMaskMinEta: " << cfg1.edgeMaskMinEta);
  ACTS_DEBUG("hitShareThreshold: " << cfg1.hitShareThreshold);
  ACTS_DEBUG("maxEndcapClusterWidth: " << cfg1.maxEndcapClusterWidth);
  ACTS_DEBUG("validateTriplets: " << cfg1.validateTriplets);
  ACTS_DEBUG("useAdaptiveCuts: " << cfg1.useAdaptiveCuts);
  ACTS_DEBUG("addTriplets: " << cfg1.addTriplets);
  ACTS_DEBUG("tauRatioCorr: " << cfg1.tauRatioCorr);
  ACTS_DEBUG("maxAbsEtaAddTriplets: " << cfg1.maxAbsEtaAddTriplets);
  ACTS_DEBUG("d0Max: " << cfg1.d0Max);
  ACTS_DEBUG("cutDPhiMax: " << cfg1.cutDPhiMax);
  ACTS_DEBUG("cutDCurvMax: " << cfg1.cutDCurvMax);
  ACTS_DEBUG("minZ0: " << cfg1.minZ0);
  ACTS_DEBUG("maxZ0: " << cfg1.maxZ0);
  ACTS_DEBUG("maxOuterRadius: " << cfg1.maxOuterRadius);
  ACTS_DEBUG("maxSeedSplitEta: " << cfg1.maxSeedSplitEta);
  ACTS_DEBUG("maxInvRadDiff: " << cfg1.maxInvRadDiff);
  ACTS_DEBUG("===== GbtsTrackFilter =====");
  const auto &cfg2 = m_cfg.trackingFilterConfig;
  ACTS_DEBUG("sigmaMS: " << cfg2.sigmaMS);
  ACTS_DEBUG("radLen: " << cfg2.radLen);
  ACTS_DEBUG("sigmaX: " << cfg2.sigmaX);
  ACTS_DEBUG("sigmaY: " << cfg2.sigmaY);
  ACTS_DEBUG("weightX: " << cfg2.weightX);
  ACTS_DEBUG("weightY: " << cfg2.weightY);
  ACTS_DEBUG("maxDChi2X: " << cfg2.maxDChi2X);
  ACTS_DEBUG("maxDChi2Y: " << cfg2.maxDChi2Y);
  ACTS_DEBUG("addHit: " << cfg2.addHit);
  ACTS_DEBUG("maxCurvature: " << cfg2.maxCurvature);
  ACTS_DEBUG("maxZ0: " << cfg2.maxZ0);
  ACTS_DEBUG("================================");
}

}  // namespace ActsExamples
