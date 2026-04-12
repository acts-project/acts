// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <iosfwd>
#include <map>
#include <memory>
#include <vector>

namespace Acts::Experimental {

/// Connection between two GBTS layers with binning information.
struct GbtsLayerConnection {
  /// @param src_ Source layer index
  /// @param dst_ Destination layer index
  GbtsLayerConnection(std::uint32_t src_, std::uint32_t dst_)
      : src(src_), dst(dst_) {};

  /// Source and destination layer indices
  std::uint32_t src{};
  /// Destination layer index
  std::uint32_t dst{};

  /// Binning table for the connection
  std::vector<std::int32_t> binTable;
};

/// Loader and container for GBTS layer connection data.
struct GbtsLayerConnectionMap {
 public:
  /// Create a GbtsLayerConnectionMap from an input stream
  /// @param inStream Input stream containing the connection data
  /// @param lrtMode Enable LRT (Large Radius Tracking) mode
  /// @return A GbtsLayerConnectionMap instance populated with the data from the stream
  static GbtsLayerConnectionMap fromStream(std::istream& inStream,
                                           bool lrtMode);
  /// Create a GbtsLayerConnectionMap from a file
  /// @param inFile Input configuration file path
  /// @param lrtMode Enable LRT (Large Radius Tracking) mode
  /// @return A GbtsLayerConnectionMap instance populated with the data from the file
  static GbtsLayerConnectionMap fromFile(std::string& inFile, bool lrtMode);

  /// Group of connections targeting a destination layer.
  struct LayerGroup {
    /// @param dst_ Destination layer key
    /// @param sources_ Vector of source connections
    LayerGroup(std::uint32_t dst_,
               const std::vector<const GbtsLayerConnection*>& sources_)
        : dst(dst_), sources(sources_) {};

    /// The target layer of the group
    std::uint32_t dst{};

    /// The source layers of the group
    std::vector<const GbtsLayerConnection*> sources;
  };

  /// Eta bin width
  float etaBinWidth{};

  /// Map of layer groups indexed by layer
  std::map<std::int32_t, std::vector<LayerGroup>> layerGroups;
  /// Map of connections indexed by layer
  std::map<std::int32_t, std::vector<std::unique_ptr<GbtsLayerConnection>>>
      connectionMap;
};

}  // namespace Acts::Experimental
