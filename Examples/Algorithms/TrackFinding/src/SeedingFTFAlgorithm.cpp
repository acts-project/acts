// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingFTFAlgorithm.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <vector>

using namespace std;

template class Acts::TrigFTF_GNN_Layer<ActsExamples::SimSpacePoint>;
template class Acts::TrigFTF_GNN_Geometry<ActsExamples::SimSpacePoint>;
template class Acts::TrigFTF_GNN_Node<ActsExamples::SimSpacePoint>;
template class Acts::TrigFTF_GNN_EtaBin<ActsExamples::SimSpacePoint>;
template struct Acts::FTF_SP<ActsExamples::SimSpacePoint>;
template class Acts::TrigFTF_GNN_DataStorage<ActsExamples::SimSpacePoint>;
template class Acts::TrigFTF_GNN_Edge<ActsExamples::SimSpacePoint>;

// constructor:
ActsExamples::SeedingFTFAlgorithm::SeedingFTFAlgorithm(
    ActsExamples::SeedingFTFAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("SeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {
  // fill config struct
  m_cfg.layerMappingFile = m_cfg.layerMappingFile;

  m_cfg.seedFilterConfig = m_cfg.seedFilterConfig.toInternalUnits();

  m_cfg.seedFinderConfig =
      m_cfg.seedFinderConfig.toInternalUnits().calculateDerivedQuantities();

  m_cfg.seedFinderOptions =
      m_cfg.seedFinderOptions.toInternalUnits().calculateDerivedQuantities(
          m_cfg.seedFinderConfig);

  for (const auto &spName : m_cfg.inputSpacePoints) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto &handle = m_inputSpacePoints.emplace_back(
        std::make_unique<ReadDataHandle<SimSpacePointContainer>>(
            this,
            "InputSpacePoints#" + std::to_string(m_inputSpacePoints.size())));
    handle->initialize(spName);
  }

  m_outputSeeds.initialize(m_cfg.outputSeeds);

  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);

  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
          Acts::SeedFilter<SimSpacePoint>(m_cfg.seedFilterConfig));

  // map
  m_cfg.ACTS_FTF_Map = Make_ACTS_FTF_Map();
  // input trig vector
  m_cfg.seedFinderConfig.m_layerGeometry = LayerNumbering();

  std::ifstream input_ifstream(
      m_cfg.seedFinderConfig.fastrack_input_file.c_str(), std::ifstream::in);

  // fastrack
  std::unique_ptr<Acts::FasTrackConnector> input_fastrack =
      std::make_unique<Acts::FasTrackConnector>(input_ifstream);

  mGNNgeo = std::make_unique<Acts::TrigFTF_GNN_Geometry<SimSpacePoint>>(
      m_cfg.seedFinderConfig.m_layerGeometry, input_fastrack);

}  // this is not FTF config type because it is a member of the algs config,
   // which is of type FTF cofig

// execute:
ActsExamples::ProcessCode ActsExamples::SeedingFTFAlgorithm::execute(
    const AlgorithmContext &ctx) const {
  std::vector<Acts::FTF_SP<SimSpacePoint>> FTF_spacePoints =
      Make_FTF_spacePoints(ctx, m_cfg.ACTS_FTF_Map);

  for (auto sp : FTF_spacePoints) {
    ACTS_DEBUG("FTF space points: "
               << " FTF_id: " << sp.FTF_ID << " z: " << sp.SP->z()
               << " r: " << sp.SP->r() << " ACTS volume:  "
               << sp.SP->sourceLinks()
                      .front()
                      .get<IndexSourceLink>()
                      .geometryId()
                      .volume()
               << "\n");
  }

  // this is now calling on a core algorithm
  Acts::SeedFinderFTF<SimSpacePoint> finder(m_cfg.seedFinderConfig, *mGNNgeo);

  // need this function as create_coords is needed for seeds
  std::function<std::pair<Acts::Vector3, Acts::Vector2>(
      const SimSpacePoint *sp)>
      create_coordinates = [](const SimSpacePoint *sp) {
        Acts::Vector3 position(sp->x(), sp->y(), sp->z());
        Acts::Vector2 variance(sp->varianceR(), sp->varianceZ());
        return std::make_pair(position, variance);
      };
  // output of function needed for seed

  finder.loadSpacePoints(FTF_spacePoints);

  Acts::RoiDescriptor internalRoi(0, -4.5, 4.5, 0, -M_PI, M_PI, 0, -150.0,
                                  150.0);

  finder.createSeeds(internalRoi, *mGNNgeo);

  // still to develop
  SimSeedContainer seeds = finder.createSeeds_old(
      m_cfg.seedFinderOptions, FTF_spacePoints, create_coordinates);

  m_outputSeeds(ctx, std::move(seeds));

  return ActsExamples::ProcessCode::SUCCESS;
}

std::map<std::pair<int, int>, std::pair<int, int>>
ActsExamples::SeedingFTFAlgorithm::Make_ACTS_FTF_Map() const {
  map<std::pair<int, int>, std::pair<int, int>> ACTS_FTF;
  std::ifstream data(
      m_cfg.layerMappingFile);  // 0 in this file refers to no FTF ID
  std::string line;
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
  // file in format ACTS_vol,ACTS_lay,ACTS_mod,FTF_id
  for (auto i : parsedCsv) {
    int ACTS_vol = stoi(i[0]);
    int ACTS_lay = stoi(i[1]);
    int ACTS_mod = stoi(i[2]);
    int FTF = stoi(i[5]);
    int eta_mod = stoi(i[6]);
    int ACTS_joint = ACTS_vol * 100 + ACTS_lay;
    ACTS_FTF.insert({{ACTS_joint, ACTS_mod}, {FTF, eta_mod}});
  }

  return ACTS_FTF;
}

std::vector<Acts::FTF_SP<ActsExamples::SimSpacePoint>>
ActsExamples::SeedingFTFAlgorithm::Make_FTF_spacePoints(
    const AlgorithmContext &ctx,
    std::map<std::pair<int, int>, std::pair<int, int>> map) const {
  // create space point vectors
  std::vector<const ActsExamples::SimSpacePoint *> spacePoints;
  std::vector<Acts::FTF_SP<ActsExamples::SimSpacePoint>> FTF_spacePoints;
  FTF_spacePoints.reserve(
      m_inputSpacePoints.size());  // not sure if this is enough

  // for loop filling space
  for (const auto &isp : m_inputSpacePoints) {
    for (const auto &spacePoint : (*isp)(ctx)) {
      // fill original space point vector
      spacePoints.push_back(&spacePoint);

      // FTF space point vector
      // loop over space points, call on map
      const auto &source_link = spacePoint.sourceLinks();
      const auto &index_source_link =
          source_link.front().get<IndexSourceLink>();

      // warning if source link empty
      if (source_link.empty()) {
        // warning in officaial acts format
        ACTS_WARNING("warning source link vector is empty");
        continue;
      }
      int ACTS_vol_id = index_source_link.geometryId().volume();
      int ACTS_lay_id = index_source_link.geometryId().layer();
      int ACTS_mod_id = index_source_link.geometryId().sensitive();

      // dont want strips or HGTD
      if (ACTS_vol_id == 2 or ACTS_vol_id == 22 or ACTS_vol_id == 23 or
          ACTS_vol_id == 24) {
        continue;
      }

      // Search for vol, lay and module=0, if doesn't esist (end) then search
      // for full thing vol*100+lay as first number in pair then 0 or mod id
      auto ACTS_joint_id = ACTS_vol_id * 100 + ACTS_lay_id;
      auto key = std::make_pair(
          ACTS_joint_id,
          0);  // here the key needs to be pair of(vol*100+lay, 0)
      auto Find = map.find(key);

      if (Find ==
          map.end()) {  // if end then make new key of (vol*100+lay, modid)
        key = std::make_pair(ACTS_joint_id, ACTS_mod_id);  // mod ID
        Find = map.find(key);
      }

      // warning if key not in map
      if (Find == map.end()) {
        ACTS_WARNING("Key not found in FTF map for volume id: "
                     << ACTS_vol_id << " and layer id: " << ACTS_lay_id);
        continue;
      }

      // now should be pixel with FTF ID:
      int FTF_id =
          Find->second
              .first;  // new map the item is a pair so want first from it

      if (FTF_id == 0) {
        ACTS_WARNING("No assigned FTF ID for key for volume id: "
                     << ACTS_vol_id << " and layer id: " << ACTS_lay_id);
      }

      // access IDs from map
      int eta_mod = Find->second.second;
      int combined_id = FTF_id * 1000 + eta_mod;

      // fill FTF vector with current sapce point and ID
      FTF_spacePoints.emplace_back(&spacePoint, FTF_id, combined_id);
    }
  }
  ACTS_VERBOSE("Space points successfully assigned FTF ID");

  return FTF_spacePoints;
}

std::vector<Acts::TrigInDetSiLayer>
ActsExamples::SeedingFTFAlgorithm::LayerNumbering() const {
  std::vector<Acts::TrigInDetSiLayer> input_vector;
  std::vector<size_t> count_vector;

  m_cfg.trackingGeometry->visitSurfaces([this, &input_vector, &count_vector](
                                            const Acts::Surface *surface) {
    Acts::GeometryIdentifier geoid = surface->geometryId();
    auto ACTS_vol_id = geoid.volume();
    auto ACTS_lay_id = geoid.layer();
    auto mod_id = geoid.sensitive();
    auto bounds_vect = surface->bounds().values();
    auto center = surface->center(Acts::GeometryContext());

    // make bounds global
    Acts::Vector3 globalFakeMom(1, 1, 1);
    Acts::Vector2 min_bound_local =
        Acts::Vector2(bounds_vect[0], bounds_vect[1]);
    Acts::Vector2 max_bound_local =
        Acts::Vector2(bounds_vect[2], bounds_vect[3]);
    Acts::Vector3 min_bound_global = surface->localToGlobal(
        Acts::GeometryContext(), min_bound_local, globalFakeMom);
    Acts::Vector3 max_bound_global = surface->localToGlobal(
        Acts::GeometryContext(), max_bound_local, globalFakeMom);

    if (min_bound_global(0) >
        max_bound_global(0)) {  // checking that not wrong way round
      min_bound_global.swap(max_bound_global);
    }

    float rc = 0.0;
    float minBound = 100000.0;
    float maxBound = -100000.0;

    // convert to FTF ID
    auto ACTS_joint_id = ACTS_vol_id * 100 + ACTS_lay_id;
    auto key =
        std::make_pair(ACTS_joint_id,
                       0);  // here the key needs to be pair of(vol*100+lay, 0)
    auto Find = m_cfg.ACTS_FTF_Map.find(key);
    int FTF_id = Find->second.first;  // new map, item is pair want first
    if (Find ==
        m_cfg.ACTS_FTF_Map
            .end()) {  // if end then make new key of (vol*100+lay, modid)
      key = std::make_pair(ACTS_joint_id, mod_id);  // mod ID
      Find = m_cfg.ACTS_FTF_Map.find(key);
      FTF_id = Find->second.first;
    }

    short barrel_ec = 0;  // a variable that says if barrrel, 0 = barrel
    int eta_mod = Find->second.second;

    // assign barrel_ec depending on FTF_layer
    if (79 < FTF_id && FTF_id < 85) {  // 80s, barrel
      barrel_ec = 0;
    } else if (89 < FTF_id && FTF_id < 99) {  // 90s positive
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

    int combined_id = FTF_id * 1000 + eta_mod;
    auto current_index =
        find_if(input_vector.begin(), input_vector.end(),
                [combined_id](auto n) { return n.m_subdet == combined_id; });
    if (current_index != input_vector.end()) {  // not end so does exist
      size_t index = std::distance(input_vector.begin(), current_index);
      input_vector[index].m_refCoord += rc;
      input_vector[index].m_minBound += minBound;
      input_vector[index].m_maxBound += maxBound;
      count_vector[index] += 1;  // increase count at the index

    } else {  // end so doesn't exists
      // make new if one with FTF ID doesn't exist:
      Acts::TrigInDetSiLayer new_FTF_ID(combined_id, barrel_ec, rc, minBound,
                                        maxBound);
      input_vector.push_back(new_FTF_ID);
      count_vector.push_back(
          1);  // so the element exists and not divinding by 0
    }

    if (m_cfg.fill_module_csv) {
      fstream fout;
      fout.open("ACTS_modules.csv",
                ios::out | ios::app);  // add to file each time
      // print to csv for each module, no repeats so dont need to make
      // map for averaging
      fout << ACTS_vol_id << ", "                                  // vol
           << ACTS_lay_id << ", "                                  // lay
           << mod_id << ", "                                       // module
           << FTF_id << ","                                        // FTF id
           << eta_mod << ","                                       // eta_mod
           << center(2) << ", "                                    // z
           << sqrt(center(0) * center(0) + center(1) * center(1))  // r
           << "\n";
    }
  });

  for (long unsigned int i = 0; i < input_vector.size(); i++) {
    input_vector[i].m_refCoord = input_vector[i].m_refCoord / count_vector[i];
  }

  return input_vector;
}
