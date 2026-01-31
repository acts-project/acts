// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <memory>
#include <vector>

namespace Acts::Experimental {

struct GbtsConnection {
  GbtsConnection(std::uint32_t s, std::uint32_t d);

  std::uint32_t m_src;
  std::uint32_t m_dst;

  std::vector<std::int32_t> m_binTable;
};

class GbtsConnector {
 public:
  struct LayerGroup {
    LayerGroup(std::uint32_t l1Key, const std::vector<const GbtsConnection*>& v)
        : m_dst(l1Key), m_sources(v) {};

    std::uint32_t m_dst{};  //< the target layer of the group

    std::vector<const GbtsConnection*>
        m_sources;  //< the source layers of the group
  };

 public:
  GbtsConnector(std::string& inFile, bool lrtMode);

  float m_etaBin{};

  std::map<std::int32_t, std::vector<struct LayerGroup>> m_layerGroups;
  std::map<std::int32_t, std::vector<std::unique_ptr<GbtsConnection>>>
      m_connMap;
};

}  // namespace Acts::Experimental
