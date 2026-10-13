// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsBinConnection.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace Acts::Experimental {

/// Builds the GBTS eta bin connection table from observed track hits.
class GbtsBinConnectionTool {
 public:
  /// Struct to hold r and z bounds for a given detector layer
  struct LayerDescription {
    /// Minimum radius
    float minR{};
    /// Maximum radius
    float maxR{};
    /// Minimum z coordinate
    float minZ{};
    /// Maximum z coordinate
    float maxZ{};
    /// layer ID
    GbtsExperimentLayerId gbtsId{};
  };

  /// Configuration for the bin connection tool
  struct Config {
    /// List of detector layers
    std::vector<LayerDescription> detectorGeometry{};

    /// Symmeterize bin connection table
    bool doSymmetrization = false;
    /// Minimum probability cut applied to bin transitions, -1 keeps them all
    float probThreshold = 1e-4f;
  };

  /// @param config Tool configuration
  /// @param geometry The geometry whose eta bins are connected, built with the
  ///                 same layers and eta bin width as the seeding it trains
  /// @param logger The Acts logger
  GbtsBinConnectionTool(const Config& config,
                        std::shared_ptr<const GbtsGeometry> geometry,
                        std::unique_ptr<const Logger> logger = getDefaultLogger(
                            "GbtsBinConnectionTool", Logging::Level::INFO));

  /// converts the bins of the hits of a track to bin transitions
  /// @param track the GBTS bins of the hits of a particle, in the order it
  ///              passed them
  void addTrack(std::span<const GbtsLayerBin> track);

  /// Creates the connection table
  /// @return the bin connections, outer to inner bin
  std::vector<GbtsBinConnection> createConnectionTable() const;

 private:
  /// returns the Acts logger
  const Logger& logger() const { return *m_logger; }

  /// gets the index to the vector of detector layers via an GBTS id
  /// @param gbtsId Gbts Id of layer
  /// @return detector layer index
  std::uint32_t getIndexByGbtsId(GbtsExperimentLayerId gbtsId) const;

  /// finds the opposide layer of a symmetrical detector with a given reference
  /// layer
  /// @param layer the detector layer
  /// @return oppsite side layers gbts id
  std::optional<GbtsExperimentLayerId> oppositeSideLayer(
      GbtsExperimentLayerId layer) const;

  /// finds the mirrored bin of a symmetrical detector: the bin of the opposite
  /// side layer with the mirrored index, as the bins run with increasing eta
  /// @param bin the bin to mirror
  /// @return the bin on the opposite side
  std::optional<GbtsLayerBin> oppositeSideBin(const GbtsLayerBin& bin) const;

  /// Config for bin connection tool
  Config m_cfg;
  /// The geometry the bins belong to
  std::shared_ptr<const GbtsGeometry> m_geometry;
  /// Acts logger
  std::unique_ptr<const Acts::Logger> m_logger;
  /// number of transitions per pair of bins, inner bin first
  std::map<std::pair<GbtsLayerBin, GbtsLayerBin>, std::uint32_t> m_binPairs{};
  /// total number of tracks used to train the table on
  std::uint32_t m_totalTracks = 0;
};

}  // namespace Acts::Experimental
