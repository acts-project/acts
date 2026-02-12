// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <vector>

namespace Acts::Experimental {

/// Connection between two GBTs layers with binning information.
struct GbtsConnection {
  /// Constructor
  /// @param s Source layer index
  /// @param d Destination layer index
  GbtsConnection(std::uint32_t s, std::uint32_t d);

  /// Source and destination layer indices
  std::uint32_t m_src;
  /// Destination layer index
  std::uint32_t m_dst;

  /// Binning table for the connection
  std::vector<std::int32_t> m_binTable;
};

/// Loader and container for GBTs layer connection data.
class GbtsConnector {
 public:
  /// Group of connections targeting a destination layer.
  struct LayerGroup {
    /// Constructor
    /// @param l1Key Destination layer key
    /// @param v Vector of source connections
    LayerGroup(std::uint32_t l1Key, const std::vector<const GbtsConnection*>& v)
        : m_dst(l1Key), m_sources(v) {};

    /// The target layer of the group
    std::uint32_t m_dst{};

    /// The source layers of the group
    std::vector<const GbtsConnection*> m_sources;
  };

 public:
  /// Constructor
  /// @param inFile Input configuration file path
  /// @param lrtMode Enable LRT (Large Radius Tracking) mode
  GbtsConnector(std::string& inFile, bool lrtMode);

  /// Eta bin size
  float m_etaBin{};

  /// Map of layer groups indexed by layer
  std::map<std::int32_t, std::vector<struct LayerGroup>> m_layerGroups;
  /// Map of connections indexed by layer
  std::map<std::int32_t, std::vector<std::unique_ptr<GbtsConnection>>>
      m_connMap;
};

}  // namespace Acts::Experimental
