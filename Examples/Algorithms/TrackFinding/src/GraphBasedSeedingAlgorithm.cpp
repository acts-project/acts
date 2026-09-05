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

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ActsExamples {

namespace {

/// What the ATLAS connector file has to say: the eta bin width its layers were
/// trained with, and the layer pairs the seeder may connect.
struct ConnectorTable {
  float etaBinWidth{};
  std::vector<Acts::Experimental::GbtsLayerConnection> connections;
};

/// Read an ATLAS GBTS connector file: a `nLinks etaBinWidth` header, then per
/// link a `index stage src dst height width nEntries` line followed by a
/// height x width bin table. Only the layer ids survive - GbtsGeometry derives
/// the bin table from the layer geometry, and the stages from the connections.
///
/// @param path Path to the connector file
/// @param stripConnections Keep the strip connections instead of the pixel ones
ConnectorTable readConnectorTable(const std::string &path,
                                  bool stripConnections) {
  std::ifstream inStream(path);
  if (!inStream) {
    throw std::runtime_error("Cannot open GBTS connector file '" + path + "'");
  }

  ConnectorTable table;

  std::uint32_t nLinks{};
  inStream >> nLinks >> table.etaBinWidth;

  // the file's own stage column, which only fixes the order the connections
  // are handed over in
  std::vector<std::pair<std::uint32_t, Acts::Experimental::GbtsLayerConnection>>
      staged;
  staged.reserve(nLinks);

  for (std::uint32_t l = 0; l < nLinks; l++) {
    std::uint32_t lIdx{};
    std::uint32_t stage{};
    std::uint32_t src{};
    std::uint32_t dst{};
    std::uint32_t height{};
    std::uint32_t width{};
    std::uint32_t nEntries{};

    inStream >> lIdx >> stage >> src >> dst >> height >> width >> nEntries;

    std::uint32_t dummy{};
    for (std::uint32_t i = 0; i < height * width; ++i) {
      inStream >> dummy;
    }

    // ATLAS ITk volume ids: 12, 13 and 14 are the strip subdetectors. The
    // table holds both technologies and only one of them is ever seeded.
    const auto isStrip = [](std::uint32_t layerId) {
      const std::uint32_t volumeId = layerId / 1000;
      return volumeId == 12 || volumeId == 13 || volumeId == 14;
    };
    if (isStrip(src) != stripConnections || isStrip(dst) != stripConnections) {
      continue;
    }

    staged.emplace_back(stage,
                        Acts::Experimental::GbtsLayerConnection{src, dst});
  }

  if (!inStream) {
    throw std::runtime_error("Malformed GBTS connector file '" + path + "'");
  }

  std::ranges::stable_sort(staged, {},
                           [](const auto &entry) { return entry.first; });

  table.connections.reserve(staged.size());
  for (const auto &[stage, connection] : staged) {
    table.connections.push_back(connection);
  }

  return table;
}

/// Read an ATLAS GBTS tau lookup table: per line a cluster width, the bulk tau
/// bounds and the near-edge ones. The width is dropped - a row is located by
/// index, one row per `tauLutBinWidth` of cluster width, never searched.
///
/// @param path Path to the lookup table file
Acts::Experimental::detail::GbtsTauLookupTable readTauLookupTable(
    const std::filesystem::path &path) {
  std::ifstream inStream(path);
  if (!inStream) {
    throw std::runtime_error("Cannot open GBTS tau lookup table '" +
                             path.string() + "'");
  }

  Acts::Experimental::detail::GbtsTauLookupTable tauLut;

  float clusterWidth{};
  Acts::Experimental::detail::GbtsTauBounds bounds;
  while (inStream >> clusterWidth >> bounds.minTau >> bounds.maxTau >>
         bounds.minTauNearEdge >> bounds.maxTauNearEdge) {
    tauLut.push_back(bounds);
  }

  if (!inStream.eof()) {
    // ended on a parse error, not on a clean EOF
    throw std::runtime_error("Malformed GBTS tau lookup table '" +
                             path.string() + "'");
  }

  return tauLut;
}

}  // namespace

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
  const ConnectorTable connectorTable = readConnectorTable(
      m_cfg.connectorInputFile, m_cfg.seedFinderConfig.useStripConnections);

  // the cluster width cuts are the only user of the tau lookup table
  if (m_cfg.seedFinderConfig.useClusterWidthCuts) {
    m_cfg.seedFinderConfig.tauLookupTable =
        readTauLookupTable(m_cfg.lutInputFile);
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

  const Acts::Experimental::GraphBasedTrackSeeder::Options options(
      m_cfg.bFieldInZ);

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

  // prepare the acts to gbts mapping file
  // 0 in this file refers to no Gbts ID
  std::ifstream data(m_cfg.layerMappingFile);
  std::string line;
  // row = physical module, column = ACTS ID components
  std::vector<std::vector<std::string>> parsedCsv;
  while (std::getline(data, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    std::vector<std::string> parsedRow;
    while (std::getline(lineStream, cell, ',')) {
      parsedRow.push_back(cell);
    }

    parsedCsv.push_back(parsedRow);
  }

  // file in format ACTS_vol,ACTS_lay,ACTS_mod,gbtsId
  for (auto i : parsedCsv) {
    const auto actsVol = static_cast<std::uint32_t>(std::stoul(i[0]));
    const auto actsLay = static_cast<std::uint32_t>(std::stoul(i[1]));
    const auto actsMod = static_cast<std::uint32_t>(std::stoul(i[2]));
    const auto gbts = static_cast<Acts::Experimental::GbtsExperimentLayerId>(
        std::stoul(i[5]));
    const auto etaMod = static_cast<std::uint32_t>(std::stoul(i[6]));
    const std::uint32_t actsJoint = actsVol * 100 + actsLay;
    const ActsIDs actsId{actsJoint, actsMod};
    const GbtsIDs gbtsId{.layerId = gbts, .etaModule = etaMod};
    actsToGbtsMap.insert({actsId, gbtsId});
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

  // dont want strips or HGTD
  if (actsVolId == 2 || actsVolId == 22 || actsVolId == 23 || actsVolId == 24) {
    return std::nullopt;
  }

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

  // warning if key not in map
  if (find == m_actsGbtsMap.end()) {
    ACTS_WARNING("Key not found in Gbts map for volume id: "
                 << actsVolId << " and layer id: " << actsLayId);
    return std::nullopt;
  }

  // now should be pixel with Gbts ID
  if (find->second.layerId == 0) {
    ACTS_WARNING("No assigned Gbts ID for key for volume id: "
                 << actsVolId << " and layer id: " << actsLayId);
  }

  return find->second.layerIndex;
}

void GraphBasedSeedingAlgorithm::resolveLayerIndices(
    const Acts::Experimental::GbtsGeometry &geometry) {
  for (auto &[actsId, gbtsId] : m_actsGbtsMap) {
    const Acts::Experimental::GbtsExperimentLayerId combinedId =
        gbtsId.layerId * 1000 + gbtsId.etaModule;

    const std::optional<Acts::Experimental::GbtsLayerIndex> index =
        geometry.layerIndex(combinedId);

    if (!index.has_value()) {
      ACTS_WARNING("No GBTS layer for combined ID: " << combinedId);
    }

    gbtsId.layerIndex = index;
  }
}

std::vector<Acts::Experimental::GbtsLayerDescription>
GraphBasedSeedingAlgorithm::layerNumbering(
    const Acts::GeometryContext &gctx) const {
  std::vector<Acts::Experimental::GbtsLayerDescription> inputVector;
  std::vector<std::size_t> countVector;

  m_cfg.trackingGeometry->visitSurfaces([this, &inputVector, &countVector,
                                         &gctx](const Acts::Surface *surface) {
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

    if (find == m_actsGbtsMap.end()) {
      ACTS_WARNING("Key not found in Gbts map for volume id: "
                   << actsVolId << ", layer id: " << actsLayId
                   << ", sensitive id: " << mod_id);
      return;  // skip this surface in the visitor
    }

    const Acts::Experimental::GbtsExperimentLayerId gbtsId =
        find->second.layerId;

    Acts::Experimental::GbtsLayerType barrelEc =
        Acts::Experimental::GbtsLayerType::Barrel;  // a variable that says if
                                                    // barrrel, 0 = barrel
    const std::uint32_t etaMod = find->second.etaModule;

    // assign barrelEc depending on Gbts_layer
    if (79 < gbtsId && gbtsId < 85) {  // 80s, barrel
      barrelEc = Acts::Experimental::GbtsLayerType::Barrel;
    } else if (89 < gbtsId && gbtsId < 99) {  // 90s positive
      barrelEc = Acts::Experimental::GbtsLayerType::Endcap;
    } else {  // 70s negative
      barrelEc = Acts::Experimental::GbtsLayerType::Endcap;
    }

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

    const Acts::Experimental::GbtsExperimentLayerId combinedId =
        gbtsId * 1000 + etaMod;

    const auto currentIndex =
        find_if(inputVector.begin(), inputVector.end(),
                [combinedId](auto n) { return n.id == combinedId; });
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
      // every layer the examples framework feeds GBTS is a pixel layer
      inputVector.push_back(Acts::Experimental::GbtsLayerDescription{
          .id = combinedId,
          .type = barrelEc,
          .technology = Acts::Experimental::GbtsLayerTechnology::Pixel,
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
           << etaMod << ","      // etaMod
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
