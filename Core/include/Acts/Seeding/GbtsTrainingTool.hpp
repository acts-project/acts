// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Acts::Experimental {

/// Layer connection training tool for
/// creating layer map for GBTS
class GbtsLayerConnectionTool {
 public:
  /// Struct to hold r and z bounds for a given detector layer
  struct LayerDescription {
    /// Constructor for filling detector layer information
    /// @param minR_ minimum radius of layer
    /// @param maxR_ maximum radius of layer
    /// @param minZ_ minimum z coordinate of layer
    /// @param maxZ_ maximum z coordinate of layer
    /// @param gbtsId_ GBTS id of layer
    LayerDescription(float minR_, float maxR_, float minZ_, float maxZ_,
                     std::int32_t gbtsId_);

    /// Minimum radius
    float minR{};
    /// Maximum radius
    float maxR{};
    /// Minimum z coordinate
    float minZ{};
    /// Maximum z coordinate
    float maxZ{};
    /// layer ID
    std::int32_t gbtsId{};
  };

  /// Configuration for the layer connection tool
  struct Config {
    /// List of detector layers
    std::vector<LayerDescription> detectorGeometry{};

    // tolerances used for assigning layer ID's

    /// Tolerance on minimum z value
    float zMinTol = 0.2340f;
    /// Tolerance on maximum z value
    float zMaxTol = 0.2340f;
    /// Tolerance on minimum radius value
    float rMinTol = 2.5337f;
    /// Tolerance on maximum radius value
    float rMaxTol = 2.5337f;

    /// Symmeterize layer connection table
    bool doSymmetrization = false;
    /// Minimum probability cut applied to layer transitions
    float probThreshold = -1;
  };

  /// Container for track layer hit information
  struct HitCoordinates {
    /// Radius value of hit
    float r{};
    /// z value of hit
    float z{};
  };

  /// pair of layer transitions
  using LayerIdPair = std::pair<std::int32_t, std::int32_t>;

  /// Hash id used for unordered sets and maps
  struct LayerIdPairHash {
    /// operator to allow the lookup of std::pair objects
    /// in unordered maps or sets
    /// @param pair Layer transition pair
    /// @return hash id
    std::size_t operator()(const LayerIdPair& pair) const noexcept {
      const auto h1 = std::hash<std::int32_t>{}(pair.first);
      const auto h2 = std::hash<std::int32_t>{}(pair.second);

      return h1 ^ (h2 << 1);
    }
  };

  /// Container of pairs of layer transitions
  using LayerIdPairs = std::unordered_set<LayerIdPair, LayerIdPairHash>;
  /// Map of layer pair transitions, quantifying the amount of times they occur
  using LayerIdPairMap =
      std::unordered_map<LayerIdPair, std::uint32_t, LayerIdPairHash>;

  /// Constructor for layer connection training tool
  /// @param config The training tool config
  /// @param logger The Acts logger
  explicit GbtsLayerConnectionTool(
      const Config& config,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("GbtsLayerConnectionTool", Logging::Level::INFO));

  /// converts layer hits to layer transitions
  /// @param track the layer hits of a particle
  void addTrack(const std::vector<HitCoordinates>& track);
  /// Creates the connection table
  /// @param outputFileLocation the location for the layer connection table
  GbtsLayerConnectionTool::LayerIdPairs createConnectionTable(
      const std::string& outputFileLocation) const;

 private:
  /// returns the Acts logger
  const Logger& logger() const { return *m_logger; }

  /// Finds the Gbts Coordinate of a given hits coordinate
  /// @param hit the coordinates of the particle hit on a layer
  std::optional<std::int32_t> findGbtsIdByCoord(
      const HitCoordinates& hit) const;

  /// gets the index to the vector of detector layers via an GBTS id
  /// @param gbtsId Gbts Id of layer
  std::uint32_t getIndexByGbtsId(std::int32_t gbtsId) const;

  /// finds the opposide layer of a symmetrical detector with a given reference
  /// layer
  /// @param layer the detecor layer
  std::optional<std::int32_t> oppositeSideLayer(std::int32_t layer) const;

  /// Config for layer connection tool
  Config m_cfg;
  /// Acts logger
  std::unique_ptr<const Acts::Logger> m_logger;
  /// map of layer transition pairs
  LayerIdPairMap m_layerPairs{};
  /// total number of tracks used to train the table on
  std::uint32_t m_totalTracks = 0;
};

}  // namespace Acts::Experimental
