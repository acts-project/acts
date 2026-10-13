// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsBinConnectionTool.hpp"

#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

namespace Acts::Experimental {

GbtsBinConnectionTool::GbtsBinConnectionTool(
    const Config& config, std::shared_ptr<const GbtsGeometry> geometry,
    std::unique_ptr<const Logger> logger)
    : m_cfg(config),
      m_geometry(std::move(geometry)),
      m_logger(std::move(logger)) {
  if (m_cfg.detectorGeometry.empty()) {
    throw std::runtime_error("File does not exist or could not be opened");
  }
  if (m_geometry == nullptr) {
    throw std::invalid_argument("Missing GBTS geometry of the bins");
  }
}

void GbtsBinConnectionTool::addTrack(std::span<const GbtsLayerBin> track) {
  if (track.size() < 2) {
    ACTS_WARNING("Track only has one measurement, skipping");
    return;
  }

  // update map with track bin transitions
  for (std::uint32_t id = 0; id + 1 < track.size(); id++) {
    const GbtsLayerBin& bin1 = track[id];
    const GbtsLayerBin& bin2 = track[id + 1];

    if (bin1.layer == bin2.layer) {
      ACTS_WARNING("Track transitions between same layer, skipping");

      continue;
    }

    m_binPairs[{bin1, bin2}] += 1;
  }

  m_totalTracks++;
}

std::vector<GbtsBinConnection> GbtsBinConnectionTool::createConnectionTable()
    const {
  if (m_totalTracks == 0) {
    throw std::runtime_error(
        "Warning: no tracks were added when creating connection table");
  }

  // obtain total transitions out of each inner bin (used as denominator of
  // probability)
  std::map<GbtsLayerBin, std::uint32_t> binTotals;
  for (const auto& [binPair, nTransitions] : m_binPairs) {
    binTotals[binPair.first] += nTransitions;
  }

  // find transitions that pass probability cut and add to temp container
  std::set<std::pair<GbtsLayerBin, GbtsLayerBin>> tempPairs;
  for (const auto& [binPair, nTransitions] : m_binPairs) {
    const float probability = static_cast<float>(nTransitions) /
                              static_cast<float>(binTotals.at(binPair.first));

    const bool passCut = (m_cfg.probThreshold == -1)
                             ? (probability != 0)
                             : (probability >= m_cfg.probThreshold);

    if (passCut) {
      tempPairs.emplace(binPair);
    }
  }

  // if symmetrizing connection table, add transitions that mirror found ones
  if (m_cfg.doSymmetrization) {
    const auto trainedPairs = tempPairs;
    for (const auto& [bin1, bin2] : trainedPairs) {
      // find mirrored bins
      const auto bin1Swapped = oppositeSideBin(bin1);
      const auto bin2Swapped = oppositeSideBin(bin2);

      if (!bin1Swapped || !bin2Swapped) {
        ACTS_WARNING("Cannot find oppisite side bin, skipping");
        continue;
      }

      // the set keeps a mirrored pair only once
      tempPairs.emplace(bin1Swapped.value(), bin2Swapped.value());
    }
  }

  // the transitions are inward -> outward, the connections outward -> inward
  std::vector<GbtsBinConnection> connections;
  connections.reserve(tempPairs.size());
  for (const auto& [inner, outer] : tempPairs) {
    connections.push_back({.src = outer, .dst = inner});
  }
  return connections;
}

std::uint32_t GbtsBinConnectionTool::getIndexByGbtsId(
    GbtsExperimentLayerId gbtsId) const {
  for (std::uint32_t idx = 0; idx < m_cfg.detectorGeometry.size(); idx++) {
    if (gbtsId == m_cfg.detectorGeometry[idx].gbtsId) {
      return idx;
    }
  }

  throw std::runtime_error("index not found for GBTS ID");
}

std::optional<GbtsExperimentLayerId> GbtsBinConnectionTool::oppositeSideLayer(
    GbtsExperimentLayerId layerId) const {
  const std::uint32_t layerIndex = getIndexByGbtsId(layerId);

  const auto& layer = m_cfg.detectorGeometry[layerIndex];

  bool switchedMinZ{};
  bool switchedMaxZ{};

  bool sameMaxR{};
  bool sameMinR{};

  for (const auto& switchedLayer : m_cfg.detectorGeometry) {
    switchedMinZ = (switchedLayer.minZ == -layer.maxZ);
    switchedMaxZ = (switchedLayer.maxZ == -layer.minZ);

    sameMaxR = (switchedLayer.maxR == layer.maxR);
    sameMinR = (switchedLayer.minR == layer.minR);

    if (switchedMinZ && switchedMaxZ && sameMaxR && sameMinR) {
      return switchedLayer.gbtsId;
    }
  }

  return std::nullopt;
}

std::optional<GbtsLayerBin> GbtsBinConnectionTool::oppositeSideBin(
    const GbtsLayerBin& bin) const {
  const auto oppositeLayer = oppositeSideLayer(bin.layer);
  if (!oppositeLayer.has_value()) {
    return std::nullopt;
  }

  const auto layerIndex = m_geometry->layerIndex(bin.layer);
  const auto oppositeIndex = m_geometry->layerIndex(oppositeLayer.value());
  if (!layerIndex.has_value() || !oppositeIndex.has_value()) {
    return std::nullopt;
  }

  // a mirrored layer has its bins in the opposite eta order
  const std::uint32_t numBins = m_geometry->layerBinning(*layerIndex).numBins;
  if (m_geometry->layerBinning(*oppositeIndex).numBins != numBins ||
      bin.bin >= numBins) {
    return std::nullopt;
  }

  return GbtsLayerBin{.layer = oppositeLayer.value(),
                      .bin = numBins - 1 - bin.bin};
}

}  // namespace Acts::Experimental
