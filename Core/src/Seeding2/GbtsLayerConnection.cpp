// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/GbtsLayerConnection.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <memory>
#include <ranges>
#include <set>
#include <unordered_map>

namespace Acts::Experimental {

GbtsLayerConnectionMap::GbtsLayerConnectionMap(std::string& inFile,
                                               bool lrtMode) {
  std::uint32_t nLinks{};

  std::ifstream input_ifstream(inFile.c_str(), std::ifstream::in);

  if (!input_ifstream.is_open()) {
    throw std::runtime_error("connection file not found");
  }

  input_ifstream >> nLinks >> etaBin;

  for (std::uint32_t l = 0; l < nLinks; l++) {
    std::uint32_t stage{};
    std::uint32_t lIdx{};
    std::uint32_t src{};
    std::uint32_t dst{};
    std::uint32_t nEntries{};
    std::uint32_t height{};
    std::uint32_t width{};

    input_ifstream >> lIdx >> stage >> src >> dst >> height >> width >>
        nEntries;

    auto pC = std::make_unique<GbtsLayerConnection>(src, dst);

    std::uint32_t dummy{};

    for (std::uint32_t i = 0; i < height; ++i) {
      for (std::uint32_t j = 0; j < width; ++j) {
        input_ifstream >> dummy;
      }
    }

    std::uint32_t srcvol_id = src / 1000;
    std::uint32_t dstvol_id = dst / 1000;

    bool srcIsStrip = (srcvol_id == 13 || srcvol_id == 12 || srcvol_id == 14);
    bool dstIsStrip = (dstvol_id == 13 || dstvol_id == 12 || dstvol_id == 14);
    if (lrtMode) {
      if (!srcIsStrip || !dstIsStrip) {
        continue;
      }
    } else {
      if (srcIsStrip || dstIsStrip) {
        continue;
      }
    }

    if (auto it = connectionMap.find(stage); it == connectionMap.end()) {
      std::vector<std::unique_ptr<GbtsLayerConnection>> v;
      v.push_back(std::move(pC));  // move the unique_ptr in
      connectionMap.emplace(stage,
                            std::move(v));  // move the vector into the map
    } else {
      it->second.push_back(std::move(pC));  // move into existing vector
    }
  }

  // re-arrange the connection stages

  std::list<const GbtsLayerConnection*> lConns;

  std::map<std::int32_t, std::vector<const GbtsLayerConnection*>> newConnMap;

  for (const auto& conn : connectionMap) {
    for (const auto& up : conn.second) {
      lConns.push_back(up.get());
    }
  }

  std::uint32_t stageCounter = 0;

  while (!lConns.empty()) {
    std::unordered_map<std::uint32_t, std::pair<std::int32_t, std::int32_t>>
        mCounter;  // layerKey, nDst, nSrc

    for (const auto& conn : lConns) {
      auto entryIt = mCounter.find(conn->dst);
      if (entryIt != mCounter.end()) {
        (*entryIt).second.first++;
      } else {
        std::uint32_t nDst = 1;
        std::uint32_t nSrc = 0;
        mCounter.insert(std::make_pair(conn->dst, std::make_pair(nDst, nSrc)));
      }

      entryIt = mCounter.find(conn->src);
      if (entryIt != mCounter.end()) {
        (*entryIt).second.second++;
      } else {
        std::uint32_t nDst = 0;
        std::uint32_t nSrc = 1;
        mCounter.insert(std::make_pair(conn->src, std::make_pair(nDst, nSrc)));
      }
    }

    // find layers with nSrc = 0

    std::set<std::uint32_t> zeroLayers;

    for (const auto& layerCounts : mCounter) {
      if (layerCounts.second.second != 0) {
        continue;
      }

      zeroLayers.insert(layerCounts.first);
    }

    // remove connections which use zeroLayer as destination

    std::vector<const GbtsLayerConnection*> theStage;

    std::list<const GbtsLayerConnection*>::iterator cIt = lConns.begin();

    while (cIt != lConns.end()) {
      if (zeroLayers.find((*cIt)->dst) !=
          zeroLayers.end()) {  // check if contains
        theStage.push_back(*cIt);
        cIt = lConns.erase(cIt);
        continue;
      }
      ++cIt;
    }
    newConnMap.insert(std::make_pair(stageCounter, theStage));
    stageCounter++;
  }

  // create layer groups

  std::uint32_t currentStage = 0;

  // the doublet making is done using "outside-in" approach hence the reverse
  // iterations

  for (const auto& it : std::views::reverse(newConnMap)) {
    const std::vector<const GbtsLayerConnection*>& vConn = it.second;

    // loop over links, extract all connections for the stage, group sources by
    // L1 (dst) index

    std::map<std::uint32_t, std::vector<const GbtsLayerConnection*>> l1ConnMap;

    for (const auto* conn : vConn) {
      std::uint32_t dst = conn->dst;

      std::map<std::uint32_t, std::vector<const GbtsLayerConnection*>>::iterator
          l1MapIt = l1ConnMap.find(dst);
      if (l1MapIt != l1ConnMap.end()) {
        (*l1MapIt).second.push_back(conn);
      } else {
        std::vector<const GbtsLayerConnection*> v = {conn};
        l1ConnMap.insert(std::make_pair(dst, v));
      }
    }

    std::vector<LayerGroup> lgv;

    lgv.reserve(l1ConnMap.size());

    for (const auto& l1Group : l1ConnMap) {
      lgv.push_back(LayerGroup(l1Group.first, l1Group.second));
    }

    layerGroups.insert(std::make_pair(currentStage, lgv));

    currentStage++;
  }

  newConnMap.clear();
}

}  // namespace Acts::Experimental
