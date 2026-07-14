// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsTrainingTool.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace Acts::Experimental {

GbtsLayerConnectionTool::LayerDescription::LayerDescription(
    float minR_, float maxR_, float minZ_, float maxZ_, std::int32_t gbtsId_)
    : minR(minR_), maxR(maxR_), minZ(minZ_), maxZ(maxZ_), gbtsId(gbtsId_) {}

GbtsLayerConnectionTool::GbtsLayerConnectionTool(
    const Config& config, std::unique_ptr<const Logger> logger)
    : m_cfg(config), m_logger(std::move(logger)) {
  if (m_cfg.detectorGeometry.empty()) {
    throw std::runtime_error("File does not exist or could not be opened");
  }

  // reserves
  const std::uint32_t pairReserve = m_cfg.detectorGeometry.size();
  m_layerPairs.reserve(pairReserve * (pairReserve - 1));

  // create map linked pairs of GBTS ids with number of transitions between
  // layers
  for (std::uint32_t i = 0; i < m_cfg.detectorGeometry.size(); i++) {
    for (std::uint32_t j = 0; j < m_cfg.detectorGeometry.size(); j++) {
      if (i == j) {
        continue;
      }

      LayerIdPair pair;
      pair.first = m_cfg.detectorGeometry[i].gbtsId;
      pair.second = m_cfg.detectorGeometry[j].gbtsId;

      // key = GBTS ids of pair, value = number of transitions
      m_layerPairs.emplace(pair, 0);
    }
  }
}

void GbtsLayerConnectionTool::addTrack(
    const std::vector<HitCoordinates>& track) {
  if (track.size() < 2) {
    ACTS_WARNING("Track only has one measurement, skipping");
    return;
  }

  // container for gbts IDs of the track
  std::vector<std::optional<std::int32_t>> layerGbtsIds{};
  layerGbtsIds.reserve(track.size());

  // find GBTS ids for all measurements in a track
  for (const auto& measurement : track) {
    const auto gbtsId = findGbtsIdByCoord(measurement);
    if (!gbtsId) {
      ACTS_WARNING("No Gbts Layer for coordinates with r: "
                   << measurement.r << " and z: " << measurement.z);
    }
    layerGbtsIds.emplace_back(gbtsId);
  }

  // update map with track layer transitions
  for (std::uint32_t id = 0; id + 1 < layerGbtsIds.size(); id++) {
    const auto& index1 = layerGbtsIds[id];
    const auto& index2 = layerGbtsIds[id + 1];

    // skip nonexistent layers ids
    if (!index1 || !index2) {
      continue;
    }

    if (index1.value() == index2.value()) {
      ACTS_WARNING("Track transitions between same layer, skipping");

      continue;
    }

    m_layerPairs[{index1.value(), index2.value()}] += 1;
  }

  m_totalTracks++;
}

GbtsLayerConnectionTool::LayerIdPairs
GbtsLayerConnectionTool::createConnectionTable(
    const std::string& outputFileLocation) const {
  if (m_totalTracks == 0) {
    throw std::runtime_error(
        "Warning: no tracks were added when creating connection table");
  }

  // define output text file
  std::ofstream outputFile(outputFileLocation);

  // obtain total incoming transitions for each src layer (used as denominator
  // of probability)
  std::vector<std::uint32_t> srcTotals;
  srcTotals.resize(m_cfg.detectorGeometry.size(), 0);

  for (const auto& [layerPair, nTransitions] : m_layerPairs) {
    std::uint32_t layerIndex = getIndexByGbtsId(layerPair.first);
    srcTotals[layerIndex] += nTransitions;
  }

  // find transitions that pass probability cut and add to temp container
  std::ofstream outputMap("map_output.txt");
  LayerIdPairs tempPairs;
  for (const auto& [layerPair, nTransitions] : m_layerPairs) {
    const std::uint32_t srcIndex = getIndexByGbtsId(layerPair.first);

    float probability{};

    if (srcTotals[srcIndex] == 0) {
      // avoid NAN error with 0/0
      probability = 0;
    } else {
      probability = static_cast<float>(nTransitions) / srcTotals[srcIndex];
    }
    outputMap << "src: " << layerPair.first << " dst: " << layerPair.second
              << " Transitions: " << nTransitions
              << " probability: " << probability << "\n";
    const bool passCut = (m_cfg.probThreshold == -1)
                             ? (probability != 0)
                             : (probability >= m_cfg.probThreshold);

    if (passCut) {
      tempPairs.emplace(layerPair.first, layerPair.second);
    }
  }

  // if symmetrizing connection table, add transitions that mirror found ones
  if (m_cfg.doSymmetrization) {
    for (const auto& layerPair : tempPairs) {
      // find swapped ids
      const auto srcSwappedId = oppositeSideLayer(layerPair.first);
      const auto dstSwappedId = oppositeSideLayer(layerPair.second);

      if (!srcSwappedId || !dstSwappedId) {
        ACTS_WARNING("Cannot find oppisite side layer, skipping");
        continue;
      }
      // search set to see if swapped pair has already been added
      const bool notAdded =
          (tempPairs.count({srcSwappedId.value(), dstSwappedId.value()}) == 0);

      // if not already added, add to output file
      if (notAdded) {
        tempPairs.emplace(srcSwappedId.value(), dstSwappedId.value());
      }
    }
  }

  return tempPairs;
}

std::optional<std::int32_t> GbtsLayerConnectionTool::findGbtsIdByCoord(
    const HitCoordinates& hit) const {
  for (const auto& layer : m_cfg.detectorGeometry) {
    const float zMin = layer.minZ - m_cfg.zMinTol;
    const float zMax = layer.maxZ + m_cfg.zMaxTol;
    const float rMin = layer.minR - m_cfg.rMinTol;
    const float rMax = layer.maxR + m_cfg.rMaxTol;

    if (zMin <= hit.z && hit.z <= zMax) {
      if (rMin <= hit.r && hit.r <= rMax) {
        return layer.gbtsId;
      }
    }
  }

  return std::nullopt;
}

std::uint32_t GbtsLayerConnectionTool::getIndexByGbtsId(
    std::int32_t gbtsId) const {
  for (std::uint32_t idx = 0; idx < m_cfg.detectorGeometry.size(); idx++) {
    if (gbtsId == m_cfg.detectorGeometry[idx].gbtsId) {
      return idx;
    }
  }

  throw std::runtime_error("index not found for GBTS ID");
}

std::optional<std::int32_t> GbtsLayerConnectionTool::oppositeSideLayer(
    std::int32_t layerId) const {
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

}  // namespace Acts::Experimental
