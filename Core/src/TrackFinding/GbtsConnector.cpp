// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// TODO: update to C++17 style
#include "Acts/TrackFinding/GbtsConnector.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <memory>
#include <set>
#include <unordered_map>

namespace Acts::Experimental {

GbtsConnection::GbtsConnection(std::uint32_t s, std::uint32_t d)
    : m_src(s), m_dst(d) {}

GbtsConnector::GbtsConnector(std::string& inFile, bool lrtMode) {
  m_connMap.clear();
  m_layerGroups.clear();

  std::uint32_t nLinks{};

  std::ifstream input_ifstream(inFile.c_str(), std::ifstream::in);

  if (!input_ifstream.is_open()) {
    throw std::runtime_error("connection file not found");
  }

  input_ifstream >> nLinks >> m_etaBin;

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

    auto pC = std::make_unique<GbtsConnection>(src, dst);

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

    auto it = m_connMap.find(stage);

    if (it == m_connMap.end()) {
      std::vector<std::unique_ptr<GbtsConnection>> v;
      v.push_back(std::move(pC));              // move the unique_ptr in
      m_connMap.emplace(stage, std::move(v));  // move the vector into the map
    } else {
      it->second.push_back(std::move(pC));  // move into existing vector
    }
  }

  // re-arrange the connection stages

  std::list<const GbtsConnection*> lConns;

  std::map<std::int32_t, std::vector<const GbtsConnection*>> newConnMap;

  for (const auto& conn : m_connMap) {
    for (const auto& up : conn.second) {
      lConns.push_back(up.get());
    }
  }

  std::uint32_t stageCounter = 0;

  while (!lConns.empty()) {
    std::unordered_map<std::uint32_t, std::pair<std::int32_t, std::int32_t>>
        mCounter;  // layerKey, nDst, nSrc

    for (const auto& conn : lConns) {
      auto entryIt = mCounter.find(conn->m_dst);
      if (entryIt != mCounter.end()) {
        (*entryIt).second.first++;
      } else {
        std::uint32_t nDst = 1;
        std::uint32_t nSrc = 0;
        mCounter.insert(
            std::make_pair(conn->m_dst, std::make_pair(nDst, nSrc)));
      }

      entryIt = mCounter.find(conn->m_src);
      if (entryIt != mCounter.end()) {
        (*entryIt).second.second++;
      } else {
        std::uint32_t nDst = 0;
        std::uint32_t nSrc = 1;
        mCounter.insert(
            std::make_pair(conn->m_src, std::make_pair(nDst, nSrc)));
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

    std::vector<const GbtsConnection*> theStage;

    std::list<const GbtsConnection*>::iterator cIt = lConns.begin();

    while (cIt != lConns.end()) {
      if (zeroLayers.find((*cIt)->m_dst) !=
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

  for (auto it = newConnMap.rbegin(); it != newConnMap.rend(); ++it) {
    const std::vector<const GbtsConnection*>& vConn = (*it).second;

    // loop over links, extract all connections for the stage, group sources by
    // L1 (dst) index

    std::map<std::uint32_t, std::vector<const GbtsConnection*>> l1ConnMap;

    for (const auto* conn : vConn) {
      std::uint32_t dst = conn->m_dst;

      std::map<std::uint32_t, std::vector<const GbtsConnection*>>::iterator
          l1MapIt = l1ConnMap.find(dst);
      if (l1MapIt != l1ConnMap.end()) {
        (*l1MapIt).second.push_back(conn);
      } else {
        std::vector<const GbtsConnection*> v = {conn};
        l1ConnMap.insert(std::make_pair(dst, v));
      }
    }

    std::vector<LayerGroup> lgv;

    lgv.reserve(l1ConnMap.size());

    for (const auto& l1Group : l1ConnMap) {
      lgv.push_back(LayerGroup(l1Group.first, l1Group.second));
    }

    m_layerGroups.insert(std::make_pair(currentStage, lgv));

    currentStage++;
  }

  newConnMap.clear();
}

}  // namespace Acts::Experimental
