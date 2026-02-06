// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
// Consider to moving to detail subdirectory
#include <map>
#include <string>
#include <vector>

namespace Acts::Experimental {

/// Connection between two GBTs layers with binning information.
struct GbtsConnection {
  /// Constructor
  /// @param s Source layer index
  /// @param d Destination layer index
  GbtsConnection(unsigned int s, unsigned int d);

  /// Source and destination layer indices
  unsigned int m_src;
  /// Destination layer index
  unsigned int m_dst;

  /// Binning table for the connection
  std::vector<int> m_binTable;
};

/// Loader and container for GBTs layer connection data.
class GbtsConnector {
 public:
  /// Group of connections targeting a destination layer.
  struct LayerGroup {
    /// Constructor
    /// @param l1Key Destination layer key
    /// @param v Vector of source connections
    LayerGroup(unsigned int l1Key, const std::vector<const GbtsConnection*>& v)
        : m_dst(l1Key), m_sources(v) {};

    /// The target layer of the group
    unsigned int m_dst;

    /// The source layers of the group
    std::vector<const GbtsConnection*> m_sources;
  };

 public:
  /// Constructor
  /// @param inFile Input configuration file path
  /// @param LRTmode Enable LRT (Large Radius Tracking) mode
  GbtsConnector(std::string& inFile, bool LRTmode);
  ~GbtsConnector();

  /// Eta bin size
  float m_etaBin{};

  /// Map of layer groups indexed by layer
  std::map<int, std::vector<struct LayerGroup> > m_layerGroups;
  /// Map of connections indexed by layer
  std::map<int, std::vector<GbtsConnection*> > m_connMap;
};

}  // namespace Acts::Experimental
