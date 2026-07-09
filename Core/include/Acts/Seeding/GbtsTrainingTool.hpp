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

class GbtsLayerConnectionTool {
 public:
  struct LayerDescription {
    LayerDescription(float minR_, float maxR_, float minZ_, float maxZ_,
                     std::int32_t gbtsId_);

    // r values
    float minR{};
    float maxR{};

    // z values
    float minZ{};
    float maxZ{};

    std::int32_t gbtsId{};
  };

  struct Config {
    std::vector<LayerDescription> detectorGeometry{};
    // tolerances used for assigning layer ID's
    float zMinTol = 0.2340f;
    float zMaxTol = 0.2340f;
    float rMinTol = 2.5337f;
    float rMaxTol = 2.5337f;

    bool doSymmetrization = false;
    bool useOldFormatting = false;

    float probThreshold = -1;
  };
  struct HitCoordinates {
    float r{};
    float z{};
  };

  using LayerIdPair = std::pair<std::int32_t, std::int32_t>;

  struct LayerIdPairHash {
    std::size_t operator()(const LayerIdPair& pair) const noexcept {
      const auto h1 = std::hash<std::int32_t>{}(pair.first);
      const auto h2 = std::hash<std::int32_t>{}(pair.second);

      return h1 ^ (h2 << 1);
    }
  };

  using LayerIdPairs = std::unordered_set<LayerIdPair, LayerIdPairHash>;
  using LayerIdPairMap =
      std::unordered_map<LayerIdPair, std::uint32_t, LayerIdPairHash>;

  explicit GbtsLayerConnectionTool(
      Config& config,
      std::unique_ptr<const Logger> logger = Acts::getDefaultLogger(
          "GbtsLayerConnectionTool", Acts::Logging::Level::INFO));

  void addTrack(const std::vector<HitCoordinates>& Track);

  void createConnectionTable(const std::string& outputFileLocations) const;

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  std::optional<std::int32_t> findGbtsIdByCoord(float r, float z) const;

  std::uint32_t getIndexByGbtsId(std::int32_t gbtsId) const;

  std::optional<std::int32_t> oppositeSideLayer(std::int32_t layer) const;

  void oldStyleFormatting(const std::string& outputFileLocations,
                          const LayerIdPairs& tempTable) const;

  Config m_cfg;

  std::unique_ptr<const Acts::Logger> m_logger =
      Acts::getDefaultLogger("Finder", Acts::Logging::Level::INFO);

  LayerIdPairMap m_layerPairs{};

  std::uint32_t m_totalTracks = 0;
};

}  // namespace Acts::Experimental
