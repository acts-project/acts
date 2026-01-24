// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/GbtsSeedingAlgorithm.hpp"

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <sstream>
#include <vector>

ActsExamples::GbtsSeedingAlgorithm::GbtsSeedingAlgorithm(
    ActsExamples::GbtsSeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("SeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {
  // initialise the spacepoint, seed and cluster handles
  m_inputSpacePoints.initialize(
      m_cfg.inputSpacePoints);  // TO DO: change bindings so it only gives a
                                // string instead of a vector
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_inputClusters.initialize(m_cfg.inputClusters);

  // parse the mapping file and turn into map
  m_cfg.ActsGbtsMap = makeActsGbtsMap();

  // create the TrigInDetSiLayers (Logical Layers),
  // as well as a map that tracks there index in m_layerGeometry
  m_layerGeometry = LayerNumbering();

  // create the connection objects
  m_connector = std::make_unique<Acts::Experimental::GbtsConnector>(
      m_cfg.seedFinderConfig.connectorInputFile,
      m_cfg.seedFinderConfig.LRTmode);

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
  m_cfg.seedFinderConfig.minPt = m_cfg.seedFinderConfig.minPt * 1000;

  m_finder = std::make_unique<Acts::Experimental::SeedFinderGbts>(
      m_cfg.seedFinderConfig, std::move(m_gbtsGeo), &m_layerGeometry,
      logger().cloneWithSuffix("GbtsFinder"));

  printSeedFinderGbtsConfig(m_cfg.seedFinderConfig);
}

// execute:
ActsExamples::ProcessCode ActsExamples::GbtsSeedingAlgorithm::execute(
    const AlgorithmContext &ctx) const {
  // take spacepoints, add variables needed for GBTS and add them to new
  // container due to how spacepoint container works, we need to keep the
  // container and the external columns we added alive this is done by using a
  // tuple of the core container and the two extra columns
  auto SpContainerComponents = MakeSpContainer(ctx, m_cfg.ActsGbtsMap);

  // used to reserve size of nodes 2D vector in core
  int max_layers = m_LayeridMap.size();

  // ROI file:Defines what region in detector we are interested in, currently
  // set to entire detector
  //  Acts::Experimental::RoiDescriptor internalRoi(0, -5, 5, 0,
  //  -std::numbers::pi, std::numbers::pi, 0, -225., 225.);
  Acts::Experimental::RoiDescriptor internalRoi(
      0, -4.5, 4.5, 0, -std::numbers::pi, std::numbers::pi, 0, -150., 150.);

  // create the seeds

  Acts::SeedContainer2 seeds =
      m_finder->createSeeds(internalRoi, SpContainerComponents, max_layers);

  // move seeds to simseedcontainer to be used down stream taking fist middle
  // and last sps currently as simseeds need to be hard types so only 3
  // spacepoint can be added but in future we should be able to have any length
  // seed
  SimSeedContainer seedContainerForStorage;
  seedContainerForStorage.reserve(seeds.size());
  for (const auto &seed : seeds) {
    auto sps = seed.spacePointIndices();
    unsigned int indices = sps.size() - 1;
    std::size_t mid = static_cast<std::size_t>(std::round(indices / 2.0));
    seedContainerForStorage.emplace_back(
        *std::get<0>(SpContainerComponents)
             .at(sps[0])  // first spacepoint
             .sourceLinks()[0]
             .get<const SimSpacePoint *>(),
        *std::get<0>(SpContainerComponents)
             .at(sps[mid])  // middle spacepoint
             .sourceLinks()[0]
             .get<const SimSpacePoint *>(),
        *std::get<0>(SpContainerComponents)
             .at(sps[indices])  // last spacepoint
             .sourceLinks()[0]
             .get<const SimSpacePoint *>());

    seedContainerForStorage.back().setVertexZ(seed.vertexZ());
    seedContainerForStorage.back().setQuality(seed.quality());
  }

  m_outputSeeds(ctx, std::move(seedContainerForStorage));

  return ActsExamples::ProcessCode::SUCCESS;
}

std::map<std::pair<int, int>, std::tuple<int, int, int>>
ActsExamples::GbtsSeedingAlgorithm::makeActsGbtsMap() const {
  std::map<std::pair<int, int>, std::tuple<int, int, int>> ActsGbts;

  // prepare the acts to gbts mapping file
  std::ifstream data(
      m_cfg.layerMappingFile);  // 0 in this file refers to no Gbts ID
  std::string line;
  std::vector<std::vector<std::string>>
      parsedCsv;  // row = physical module, column = ACTS ID components
  while (std::getline(data, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    std::vector<std::string> parsedRow;
    while (std::getline(lineStream, cell, ',')) {
      parsedRow.push_back(cell);
    }

    parsedCsv.push_back(parsedRow);
  }

  // file in format ACTS_vol,ACTS_lay,ACTS_mod,Gbts_id
  for (auto i : parsedCsv) {
    int ACTS_vol = stoi(i[0]);
    int ACTS_lay = stoi(i[1]);
    int ACTS_mod = stoi(i[2]);
    int Gbts = stoi(i[5]);
    int eta_mod = stoi(i[6]);
    int ACTS_joint = ACTS_vol * 100 + ACTS_lay;
    ActsGbts.insert({{ACTS_joint, ACTS_mod}, {Gbts, eta_mod, 0}});
  }

  return ActsGbts;
}

Acts::Experimental::SPContainerComponentsType
ActsExamples::GbtsSeedingAlgorithm::MakeSpContainer(
    const AlgorithmContext &ctx,
    std::map<std::pair<int, int>, std::tuple<int, int, int>> map) const {
  // new seeding container test
  // initialise input spacepoints from handle and define new container
  const SimSpacePointContainer &spacePoints = m_inputSpacePoints(ctx);
  Acts::SpacePointContainer2 coreSpacePoints(

      Acts::SpacePointColumns::SourceLinks | Acts::SpacePointColumns::X |
      Acts::SpacePointColumns::Y | Acts::SpacePointColumns::Z |
      Acts::SpacePointColumns::R | Acts::SpacePointColumns::Phi);

  // add new column for layer ID and clusterwidth
  auto LayerColoumn = coreSpacePoints.createColumn<int>("LayerID");
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

    int ACTS_vol_id = indexSourceLink.geometryId().volume();
    int ACTS_lay_id = indexSourceLink.geometryId().layer();
    int ACTS_mod_id = indexSourceLink.geometryId().sensitive();

    // dont want strips or HGTD
    if (ACTS_vol_id == 2 || ACTS_vol_id == 22 || ACTS_vol_id == 23 ||
        ACTS_vol_id == 24) {
      continue;
    }

    // Search for vol, lay and module=0, if doesn't esist (end) then search
    // for full thing vol*100+lay as first number in pair then 0 or mod id
    auto ACTS_joint_id = ACTS_vol_id * 100 + ACTS_lay_id;
    auto key =
        std::make_pair(ACTS_joint_id,
                       0);  // here the key needs to be pair of(vol*100+lay, 0)
    auto Find = map.find(key);

    if (Find ==
        map.end()) {  // if end then make new key of (vol*100+lay, modid)
      key = std::make_pair(ACTS_joint_id, ACTS_mod_id);  // mod ID
      Find = map.find(key);
    }

    // warning if key not in map
    if (Find == map.end()) {
      ACTS_WARNING("Key not found in Gbts map for volume id: "
                   << ACTS_vol_id << " and layer id: " << ACTS_lay_id);
      continue;
    }

    // now should be pixel with Gbts ID:
    int Gbts_id = std::get<0>(
        Find->second);  // new map the item is a pair so want first from it

    if (Gbts_id == 0) {
      ACTS_WARNING("No assigned Gbts ID for key for volume id: "
                   << ACTS_vol_id << " and layer id: " << ACTS_lay_id);
    }

    // access IDs from map

    auto newSp = coreSpacePoints.createSpacePoint();

    newSp.assignSourceLinks(
        std::array<Acts::SourceLink, 1>{Acts::SourceLink(&spacePoint)});
    newSp.x() = spacePoint.x();
    newSp.y() = spacePoint.y();
    newSp.z() = spacePoint.z();
    newSp.r() = spacePoint.r();
    newSp.phi() = std::atan2(spacePoint.y(), spacePoint.x());
    newSp.extra(LayerColoumn) = std::get<2>(Find->second);
    // false input as this is not available in examples
    newSp.extra(ClusterWidthColoumn) = 0;
    newSp.extra(LocalPositionColoumn) = 0;
  }

  ACTS_VERBOSE("Space point collection successfully assigned LayerID's");

  return std::make_tuple(std::move(coreSpacePoints), LayerColoumn.asConst(),
                         ClusterWidthColoumn.asConst(),
                         LocalPositionColoumn.asConst());
}

std::vector<Acts::Experimental::TrigInDetSiLayer>
ActsExamples::GbtsSeedingAlgorithm::LayerNumbering() const {
  std::vector<Acts::Experimental::TrigInDetSiLayer> input_vector{};
  std::vector<std::size_t> count_vector{};

  m_cfg.trackingGeometry->visitSurfaces([this, &input_vector, &count_vector](
                                            const Acts::Surface *surface) {
    Acts::GeometryIdentifier geoid = surface->geometryId();
    auto ACTS_vol_id = geoid.volume();
    auto ACTS_lay_id = geoid.layer();
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

    if (min_bound_global(0) >
        max_bound_global(0)) {  // checking that not wrong way round
      min_bound_global.swap(max_bound_global);
    }

    float rc = 0.0;
    float minBound = 100000.0;
    float maxBound = -100000.0;

    // convert to Gbts ID
    auto ACTS_joint_id = ACTS_vol_id * 100 + ACTS_lay_id;
    auto key =
        std::make_pair(ACTS_joint_id,
                       0);  // here the key needs to be pair of(vol*100+lay, 0)
    auto Find = m_cfg.ActsGbtsMap.find(key);
    int Gbts_id = 0;  // initialise first to avoid FLTUND later
    Gbts_id = std::get<0>(Find->second);  // new map, item is pair want first
    if (Find ==
        m_cfg.ActsGbtsMap
            .end()) {  // if end then make new key of (vol*100+lay, modid)
      key = std::make_pair(ACTS_joint_id, mod_id);  // mod ID
      Find = m_cfg.ActsGbtsMap.find(key);
      Gbts_id = std::get<0>(Find->second);
    }

    short barrel_ec = 0;  // a variable that says if barrrel, 0 = barrel
    int eta_mod = std::get<1>(Find->second);

    // assign barrel_ec depending on Gbts_layer
    if (79 < Gbts_id && Gbts_id < 85) {  // 80s, barrel
      barrel_ec = 0;
    } else if (89 < Gbts_id && Gbts_id < 99) {  // 90s positive
      barrel_ec = 2;
    } else {  // 70s negative
      barrel_ec = -2;
    }

    if (barrel_ec == 0) {
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

    int combined_id = Gbts_id * 1000 + eta_mod;

    auto current_index =
        find_if(input_vector.begin(), input_vector.end(),
                [combined_id](auto n) { return n.m_subdet == combined_id; });
    if (current_index != input_vector.end()) {  // not end so does exist
      std::size_t index = std::distance(input_vector.begin(), current_index);
      input_vector[index].m_refCoord += rc;
      input_vector[index].m_minBound += minBound;
      input_vector[index].m_maxBound += maxBound;
      count_vector[index] += 1;  // increase count at the index

    } else {  // end so doesn't exists
      // make new if one with Gbts ID doesn't exist:
      Acts::Experimental::TrigInDetSiLayer new_Gbts_ID(combined_id, barrel_ec,
                                                       rc, minBound, maxBound);
      input_vector.push_back(new_Gbts_ID);
      count_vector.push_back(
          1);  // so the element exists and not divinding by 0

      // tracking the index of each TrigInDetSiLayer as there added to the
      // vector
      int LayerID =
          count_vector.size() -
          1;  // so layer ID refers to actual index and not size of vector
      std::get<2>(Find->second) = LayerID;
      m_LayeridMap.insert({combined_id, LayerID});
    }
    // look up for every combined ID to see if it has a layer
    auto FindLayer = m_LayeridMap.find(combined_id);
    if (FindLayer == m_LayeridMap.end()) {
      ACTS_WARNING("No assigned Layer ID for combined ID: " << combined_id);
    } else {
      std::get<2>(Find->second) = FindLayer->second;
    }

    if (m_cfg.fill_module_csv) {
      std::fstream fout;
      fout.open("ACTS_modules.csv",
                std::ios::out | std::ios::app);  // add to file each time
      // print to csv for each module, no repeats so dont need to make
      // map for averaging
      fout << ACTS_vol_id << ", "                                  // vol
           << ACTS_lay_id << ", "                                  // lay
           << mod_id << ", "                                       // module
           << Gbts_id << ","                                       // Gbts id
           << eta_mod << ","                                       // eta_mod
           << center(2) << ", "                                    // z
           << sqrt(center(0) * center(0) + center(1) * center(1))  // r
           << "\n";
    }
  });

  for (std::size_t i = 0; i < input_vector.size(); i++) {
    input_vector[i].m_refCoord = input_vector[i].m_refCoord / count_vector[i];
  }

  return input_vector;
}

void ActsExamples::GbtsSeedingAlgorithm::printSeedFinderGbtsConfig(
    const Acts::Experimental::SeedFinderGbtsConfig &cfg) {
  ACTS_DEBUG("===== SeedFinderGbtsConfig =====");

  ACTS_DEBUG("BeamSpotCorrection: " << cfg.BeamSpotCorrection
                                    << " (default: false)");
  ACTS_DEBUG("connectorInputFile: " << cfg.connectorInputFile
                                    << " (default: empty string)");
  ACTS_DEBUG("lutInputFile: " << cfg.lutInputFile
                              << " (default: empty string)");
  ACTS_DEBUG("LRTmode: " << cfg.LRTmode << " (default: false)");
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
