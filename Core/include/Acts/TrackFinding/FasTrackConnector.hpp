// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
// Consider to moving to detail subdirectory
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace Acts {

struct FasTrackConnection {
 public:
  FasTrackConnection(unsigned int s, unsigned int d);

  unsigned int m_src, m_dst;
  std::vector<int> m_binTable;
};

class FasTrackConnector {
 public:
  struct LayerGroup {
    LayerGroup(unsigned int l1Key,
               const std::vector<const Acts::FasTrackConnection *> &v)
        : m_dst(l1Key), m_sources(v) {}

    unsigned int m_dst;  // the target layer of the group
    std::vector<const Acts::FasTrackConnection *>
        m_sources;  // the source layers of the group
  };

  FasTrackConnector(std::ifstream &inFile);

  ~FasTrackConnector();

  float m_etaBin{};

  std::map<int, std::vector<struct LayerGroup>> m_layerGroups;
  std::map<int, std::vector<Acts::FasTrackConnection *>> m_connMap;
  // TODO: change to std::map<int, std::vector<Acts::FasTrackConnection> >
  // m_connMap; or   std::map<int,
  // std::vector<std::unique_ptr<Acts::FasTrackConnection>> > m_connMap;
};

}  // namespace Acts
