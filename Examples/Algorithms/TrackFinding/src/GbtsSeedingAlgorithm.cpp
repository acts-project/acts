// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/GbtsSeedingAlgorithm.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <sstream>
#include <vector>

namespace ActsExamples {

GbtsSeedingAlgorithm::GbtsSeedingAlgorithm(GbtsSeedingAlgorithm::Config cfg,
                                           Acts::Logging::Level lvl)
    : IAlgorithm("SeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {
  // initialise the spacepoint, seed and cluster handles
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_inputClusters.initialize(m_cfg.inputClusters);

  // parse the mapping file and turn into map
  m_cfg.actsGbtsMap = makeActsGbtsMap();

  // create the TrigInDetSiLayers (Logical Layers),
  // as well as a map that tracks there index in m_layerGeometry
  m_layerGeometry = layerNumbering();

  // create the connection objects
  m_connector = std::make_unique<Acts::Experimental::GbtsConnector>(
      m_cfg.seedFinderConfig.connectorInputFile,
      m_cfg.seedFinderConfig.lrtMode);

  // option that allows for adding custom eta binning (default is at 0.2)
  if (m_cfg.seedFinderConfig.etaBinOverride != 0.0f) {
    m_connector->m_etaBin = m_cfg.seedFinderConfig.etaBinOverride;
  }

  // initialise the object that holds all the geometry information needed for
  // the algorithm
  m_gbtsGeo = std::make_unique<Acts::Experimental::GbtsGeometry>(
      m_layerGeometry, m_connector);

  // manually convert min Pt as no conversion available in ACTS Examples
  // (currently inputs as 0.9 GeV but need 900 MeV)

  m_finder = std::make_unique<Acts::Experimental::SeedFinderGbts>(
      m_cfg.seedFinderConfig, std::move(m_gbtsGeo), &m_layerGeometry,
      logger().cloneWithSuffix("GbtsFinder"));

  printSeedFinderGbtsConfig(m_cfg.seedFinderConfig);
}

ProcessCode GbtsSeedingAlgorithm::execute(const AlgorithmContext &ctx) const {
  // initialise input spacepoints from handle and define new container
  const SpacePointContainer &spacePoints = m_inputSpacePoints(ctx);

  // take spacepoints, add variables needed for GBTS and add them to new
  // container due to how spacepoint container works, we need to keep the
  // container and the external columns we added alive this is done by using a
  // tuple of the core container and the two extra columns
  auto spContainerComponents = makeSpContainer(spacePoints, m_cfg.actsGbtsMap);
  const auto &coreSpacePoints =
      std::get<Acts::SpacePointContainer2>(spContainerComponents);

  // used to reserve size of nodes 2D vector in core
  std::uint32_t maxLayers = m_LayeridMap.size();

  // ROI file:Defines what region in detector we are interested in, currently
  // set to entire detector
  Acts::Experimental::RoiDescriptor internalRoi(
      0, -4.5, 4.5, 0, -std::numbers::pi, std::numbers::pi, 0, -150., 150.);

  // create the seeds
  SeedContainer seeds =
      m_finder->createSeeds(internalRoi, spContainerComponents, maxLayers);

  // update seed space point indices to original space point container
  for (auto seed : seeds) {
    for (auto &spIndex : seed.spacePointIndices()) {
      spIndex =
          coreSpacePoints.at(spIndex).sourceLinks()[0].get<SpacePointIndex>();
    }
  }

  m_outputSeeds(ctx, std::move(seeds));

  return ProcessCode::SUCCESS;
}

std::map<GbtsSeedingAlgorithm::ActsIDs, GbtsSeedingAlgorithm::GbtsIDs>
GbtsSeedingAlgorithm::makeActsGbtsMap() const {
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
    std::uint32_t actsVol = stoi(i[0]);
    std::uint32_t actsLay = stoi(i[1]);
    std::uint32_t actsMod = stoi(i[2]);
    std::uint32_t gbts = stoi(i[5]);
    std::uint32_t etaMod = stoi(i[6]);
    std::uint32_t actsJoint = actsVol * 100 + actsLay;
    ActsIDs actsId{static_cast<std::uint64_t>(actsJoint),
                   static_cast<std::uint64_t>(actsMod)};
    GbtsIDs gbtsId{static_cast<std::uint32_t>(gbts),
                   static_cast<std::uint32_t>(etaMod), 0};
    actsToGbtsMap.insert({{actsId}, {gbtsId}});
  }

  return actsToGbtsMap;
}

Acts::Experimental::SPContainerComponentsType
GbtsSeedingAlgorithm::makeSpContainer(const SpacePointContainer &spacePoints,
                                      std::map<ActsIDs, GbtsIDs> map) const {
  Acts::SpacePointContainer2 coreSpacePoints(
      Acts::SpacePointColumns::SourceLinks | Acts::SpacePointColumns::X |
      Acts::SpacePointColumns::Y | Acts::SpacePointColumns::Z |
      Acts::SpacePointColumns::R | Acts::SpacePointColumns::Phi);

  // add new column for layer ID and clusterwidth
  auto LayerColoumn = coreSpacePoints.createColumn<std::uint32_t>("layerID");
  auto ClusterWidthColoumn =
      coreSpacePoints.createColumn<float>("Cluster_Width");
  auto LocalPositionColoumn =
      coreSpacePoints.createColumn<float>("LocalPositionY");
  coreSpacePoints.reserve(spacePoints.size());

  // for loop filling space

  for (const auto &spacePoint : spacePoints) {
    // Gbts space point vector
    // loop over space points, call on map
    const auto &sourceLink = spacePoint.sourceLinks();

    // warning if source link empty
    if (sourceLink.empty()) {
      // warning in officaial acts format
      ACTS_WARNING("warning source link vector is empty");
      continue;
    }

    const auto &indexSourceLink = sourceLink.front().get<IndexSourceLink>();

    std::uint32_t actsVolId = indexSourceLink.geometryId().volume();
    std::uint32_t actsLayId = indexSourceLink.geometryId().layer();
    std::uint32_t actsModId = indexSourceLink.geometryId().sensitive();

    // dont want strips or HGTD
    if (actsVolId == 2 || actsVolId == 22 || actsVolId == 23 ||
        actsVolId == 24) {
      continue;
    }

    // Search for vol, lay and module=0, if doesn't esist (end) then search
    // for full thing vol*100+lay as first number in pair then 0 or mod id
    auto actsJointId = actsVolId * 100 + actsLayId;

    // here the key needs to be pair of(vol*100+lay, 0)
    ActsIDs key{static_cast<std::uint64_t>(actsJointId), 0};
    auto find = map.find(key);

    // if end then make new key of (vol*100+lay, modid)
    if (find == map.end()) {
      key = ActsIDs{static_cast<std::uint64_t>(actsJointId),
                    static_cast<std::uint64_t>(actsModId)};  // mod ID
      find = map.find(key);
    }

    // warning if key not in map
    if (find == map.end()) {
      ACTS_WARNING("Key not found in Gbts map for volume id: "
                   << actsVolId << " and layer id: " << actsLayId);
      continue;
    }

    // now should be pixel with Gbts ID:
    // new map the item is a pair so want first from it
    std::uint32_t gbtsId = std::get<0>(find->second);

    if (gbtsId == 0) {
      ACTS_WARNING("No assigned Gbts ID for key for volume id: "
                   << actsVolId << " and layer id: " << actsLayId);
    }

    // access IDs from map

    auto newSp = coreSpacePoints.createSpacePoint();

    newSp.assignSourceLinks(
        std::array<Acts::SourceLink, 1>{Acts::SourceLink(spacePoint.index())});
    newSp.x() = spacePoint.x();
    newSp.y() = spacePoint.y();
    newSp.z() = spacePoint.z();
    newSp.r() = spacePoint.r();
    newSp.phi() = std::atan2(spacePoint.y(), spacePoint.x());
    newSp.extra(LayerColoumn) = std::get<2>(find->second);
    // false input as this is not available in examples
    newSp.extra(ClusterWidthColoumn) = 0;
    newSp.extra(LocalPositionColoumn) = 0;
  }

  ACTS_VERBOSE("Space point collection successfully assigned layerID's");

  return std::make_tuple(std::move(coreSpacePoints), LayerColoumn.asConst(),
                         ClusterWidthColoumn.asConst(),
                         LocalPositionColoumn.asConst());
}

std::vector<Acts::Experimental::TrigInDetSiLayer>
GbtsSeedingAlgorithm::layerNumbering() const {
  std::vector<Acts::Experimental::TrigInDetSiLayer> inputVector{};
  std::vector<std::size_t> countVector{};

  m_cfg.trackingGeometry->visitSurfaces([this, &inputVector, &countVector](
                                            const Acts::Surface *surface) {
    Acts::GeometryIdentifier geoid = surface->geometryId();
    auto actsVolId = geoid.volume();
    auto actsLayId = geoid.layer();
    auto mod_id = geoid.sensitive();
    auto bounds_vect = surface->bounds().values();
    auto center =
        surface->center(Acts::GeometryContext::dangerouslyDefaultConstruct());

    // make bounds global
    Acts::Vector3 globalFakeMom(1, 1, 1);
    Acts::Vector2 min_bound_local =
        Acts::Vector2(bounds_vect[0], bounds_vect[1]);
    Acts::Vector2 max_bound_local =
        Acts::Vector2(bounds_vect[2], bounds_vect[3]);
    Acts::Vector3 min_bound_global = surface->localToGlobal(
        Acts::GeometryContext::dangerouslyDefaultConstruct(), min_bound_local,
        globalFakeMom);
    Acts::Vector3 max_bound_global = surface->localToGlobal(
        Acts::GeometryContext::dangerouslyDefaultConstruct(), max_bound_local,
        globalFakeMom);

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
    auto find = m_cfg.actsGbtsMap.find(key);
    // initialise first to avoid FLTUND later
    std::uint32_t gbtsId = 0;
    // new map, item is pair want first
    gbtsId = std::get<0>(find->second);
    // if end then make new key of (vol*100+lay, modid)
    if (find == m_cfg.actsGbtsMap.end()) {
      key = ActsIDs{actsJointId, mod_id};  // mod ID
      find = m_cfg.actsGbtsMap.find(key);
      gbtsId = std::get<0>(find->second);
    }

    std::int16_t barrelEc = 0;  // a variable that says if barrrel, 0 = barrel
    std::uint32_t etaMod = std::get<1>(find->second);

    // assign barrelEc depending on Gbts_layer
    if (79 < gbtsId && gbtsId < 85) {  // 80s, barrel
      barrelEc = 0;
    } else if (89 < gbtsId && gbtsId < 99) {  // 90s positive
      barrelEc = 2;
    } else {  // 70s negative
      barrelEc = -2;
    }

    if (barrelEc == 0) {
      rc = sqrt(center(0) * center(0) +
                center(1) * center(1));  // barrel center in r
      // bounds of z
      if (min_bound_global(2) < minBound) {
        minBound = min_bound_global(2);
      }
      if (max_bound_global(2) > maxBound) {
        maxBound = max_bound_global(2);
      }
    } else {
      rc = center(2);  // not barrel center in Z
      // bounds of r
      float min = sqrt(min_bound_global(0) * min_bound_global(0) +
                       min_bound_global(1) * min_bound_global(1));
      float max = sqrt(max_bound_global(0) * max_bound_global(0) +
                       max_bound_global(1) * max_bound_global(1));
      if (min < minBound) {
        minBound = min;
      }
      if (max > maxBound) {
        maxBound = max;
      }
    }

    std::int32_t combined_id = gbtsId * 1000 + etaMod;

    auto current_index =
        find_if(inputVector.begin(), inputVector.end(),
                [combined_id](auto n) { return n.m_subdet == combined_id; });
    if (current_index != inputVector.end()) {  // not end so does exist
      std::size_t index = std::distance(inputVector.begin(), current_index);
      inputVector[index].m_refCoord += rc;
      inputVector[index].m_minBound += minBound;
      inputVector[index].m_maxBound += maxBound;
      countVector[index] += 1;  // increase count at the index

    } else {  // end so doesn't exists
      // make new if one with Gbts ID doesn't exist:
      Acts::Experimental::TrigInDetSiLayer new_Gbts_ID(combined_id, barrelEc,
                                                       rc, minBound, maxBound);
      inputVector.push_back(new_Gbts_ID);
      // so the element exists and not divinding by 0
      countVector.push_back(1);

      // tracking the index of each TrigInDetSiLayer as there added

      // so layer ID refers to actual index and not size of vector
      std::uint32_t layerID = countVector.size() - 1;
      std::get<2>(find->second) = layerID;
      m_LayeridMap.insert({combined_id, layerID});
    }
    // look up for every combined ID to see if it has a layer
    auto FindLayer = m_LayeridMap.find(combined_id);
    if (FindLayer == m_LayeridMap.end()) {
      ACTS_WARNING("No assigned Layer ID for combined ID: " << combined_id);
    } else {
      std::get<2>(find->second) = FindLayer->second;
    }

    // add to file each time,
    // print to csv for each module, no repeats so dont need to make
    // map for averaging
    if (m_cfg.fill_module_csv) {
      std::fstream fout;
      fout.open("ACTS_modules.csv", std::ios::out | std::ios::app);
      fout << actsVolId << ", "                                    // vol
           << actsLayId << ", "                                    // lay
           << mod_id << ", "                                       // module
           << gbtsId << ","                                        // Gbts id
           << etaMod << ","                                        // etaMod
           << center(2) << ", "                                    // z
           << sqrt(center(0) * center(0) + center(1) * center(1))  // r
           << "\n";
    }
  });

  for (std::size_t i = 0; i < inputVector.size(); i++) {
    inputVector[i].m_refCoord = inputVector[i].m_refCoord / countVector[i];
  }

  return inputVector;
}

void GbtsSeedingAlgorithm::printSeedFinderGbtsConfig(
    const Acts::Experimental::SeedFinderGbtsConfig &cfg) {
  ACTS_DEBUG("===== SeedFinderGbtsConfig =====");

  ACTS_DEBUG("BeamSpotCorrection: " << cfg.BeamSpotCorrection
                                    << " (default: false)");
  ACTS_DEBUG("connectorInputFile: " << cfg.connectorInputFile
                                    << " (default: empty string)");
  ACTS_DEBUG("lutInputFile: " << cfg.lutInputFile
                              << " (default: empty string)");
  ACTS_DEBUG("lrtMode: " << cfg.lrtMode << " (default: false)");
  ACTS_DEBUG("useML: " << cfg.useML << " (default: false)");
  ACTS_DEBUG("matchBeforeCreate: " << cfg.matchBeforeCreate
                                   << " (default: false)");
  ACTS_DEBUG("useOldTunings: " << cfg.useOldTunings << " (default: false)");
  ACTS_DEBUG("tau_ratio_cut: " << cfg.tau_ratio_cut << " (default: 0.007)");
  ACTS_DEBUG("tau_ratio_precut: " << cfg.tau_ratio_precut
                                  << " (default: 0.009f)");
  ACTS_DEBUG("etaBinOverride: " << cfg.etaBinOverride << " (default: 0.0)");
  ACTS_DEBUG("nMaxPhiSlice: " << cfg.nMaxPhiSlice << " (default: 53)");
  ACTS_DEBUG("minPt: " << cfg.minPt
                       << " (default: 1.0 * Acts::UnitConstants::GeV)");
  ACTS_DEBUG("phiSliceWidth: " << cfg.phiSliceWidth << " (default: derived)");
  ACTS_DEBUG("ptCoeff: " << cfg.ptCoeff
                         << " (default: 0.29997 * 1.9972 / 2.0)");
  ACTS_DEBUG("useEtaBinning: " << cfg.useEtaBinning << " (default: true)");
  ACTS_DEBUG("doubletFilterRZ: " << cfg.doubletFilterRZ << " (default: true)");
  ACTS_DEBUG("nMaxEdges: " << cfg.nMaxEdges << " (default: 2000000)");
  ACTS_DEBUG("minDeltaRadius: " << cfg.minDeltaRadius << " (default: 2.0)");
  ACTS_DEBUG("sigmaMS: " << cfg.sigmaMS << " (default: 0.016)");
  ACTS_DEBUG("radLen: " << cfg.radLen << " (default: 0.025)");
  ACTS_DEBUG("sigma_x: " << cfg.sigma_x << " (default: 0.08)");
  ACTS_DEBUG("sigma_y: " << cfg.sigma_y << " (default: 0.25)");
  ACTS_DEBUG("weight_x: " << cfg.weight_x << " (default: 0.5)");
  ACTS_DEBUG("weight_y: " << cfg.weight_y << " (default: 0.5)");
  ACTS_DEBUG("maxDChi2_x: " << cfg.maxDChi2_x << " (default: 5.0)");
  ACTS_DEBUG("maxDChi2_y: " << cfg.maxDChi2_y << " (default: 6.0)");
  ACTS_DEBUG("add_hit: " << cfg.add_hit << " (default: 14.0)");
  ACTS_DEBUG("max_curvature: " << cfg.max_curvature << " (default: 1e-3f)");
  ACTS_DEBUG("max_z0: " << cfg.max_z0 << " (default: 170.0)");
  ACTS_DEBUG("edge_mask_min_eta: " << cfg.edge_mask_min_eta
                                   << " (default: 1.5)");
  ACTS_DEBUG("hit_share_threshold: " << cfg.hit_share_threshold
                                     << " (default: 0.49)");
  ACTS_DEBUG("max_endcap_clusterwidth: " << cfg.max_endcap_clusterwidth
                                         << " (default: 0.35)");

  ACTS_DEBUG("================================");
}

}  // namespace ActsExamples
