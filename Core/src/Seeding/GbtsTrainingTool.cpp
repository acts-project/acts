// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/GbtsTrainingTool.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

namespace Acts::Experimental {

GbtsLayerConnectionTool::LayerDescription::LayerDescription(
    float minR_, float maxR_, float minZ_, float maxZ_, std::int32_t gbtsId_)
    : minR(minR_), maxR(maxR_), minZ(minZ_), maxZ(maxZ_), gbtsId(gbtsId_) {}

GbtsLayerConnectionTool::GbtsLayerConnectionTool(
    Config config, const std::string& geometryInformation,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(config), m_logger(std::move(logger)) {
  std::ifstream inStream(geometryInformation.c_str());

  if (!inStream) {
    throw std::runtime_error("File does not exist or could not be opened");
  }

  // define how many lines there are for reserving
  std::uint32_t lines{};
  std::string line{};
  while (std::getline(inStream, line)) {
    lines++;
  }
  inStream.clear();
  inStream.seekg(0);

  // reserves
  m_detectorGeometry.reserve(lines);
  m_layerPairs.reserve(lines * (lines - 1));

  // create geometry objects
  float minR{};
  float maxR{};

  float minZ{};
  float maxZ{};

  std::int32_t gbtsId{};

  for (std::uint32_t l = 0; l < lines; l++) {
    inStream >> minR >> maxR >> minZ >> maxZ >> gbtsId;

    m_detectorGeometry.emplace_back(minR, maxR, minZ, maxZ, gbtsId);
  }

  // create map linked pairs of GBTS ids with number of transitions between
  // layers
  for (std::uint32_t i = 0; i < m_detectorGeometry.size(); i++) {
    for (std::uint32_t j = 0; j < m_detectorGeometry.size(); j++) {
      if (i == j) {
        continue;
      }

      LayerIdPair pair;
      pair.first = m_detectorGeometry[i].gbtsId;
      pair.second = m_detectorGeometry[j].gbtsId;

      // key = GBTS ids of pair, value = number of transitions
      m_layerPairs.emplace(pair, 0);
    }
  }
}

void GbtsLayerConnectionTool::addTrack(
    const std::vector<TrackCoordinates>& track) {
  if (track.size() < 2) {
    ACTS_WARNING("Track only has one measurement, skipping");
    return;
  }

  // check to see if tracks have backward transitions
  for (std::uint32_t idx = 0; idx + 1 < track.size(); idx++) {
    const float hit1Magnitude = std::hypot(track[idx].r, track[idx].z);
    const float hit2Magnitude = std::hypot(track[idx + 1].r, track[idx + 1].z);

    if (hit1Magnitude > hit2Magnitude) {
      ACTS_WARNING("Track travels inwards, skipping");
      return;
    }
  }

  // container for gbts IDs of the track
  std::vector<std::optional<std::int32_t>> layerGbtsIds{};
  layerGbtsIds.reserve(track.size());

  // find GBTS ids for all measurements in a track
  for (const auto& measurement : track) {
    const float r = measurement.r;
    const float z = measurement.z;

    const auto gbtsId = findGbtsIdByCoord(r, z);
    if (!gbtsId) {
      ACTS_WARNING("No Gbts Layer for coordinates with r: " << r
                                                            << " and z: " << z);
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

void GbtsLayerConnectionTool::createConnectionTable(
    const std::string& outputFileLocation) const {
  if (m_totalTracks == 0) {
    throw std::runtime_error(
        "Warning: no tracks were added when creating connection table");
  }

  // define output text file
  std::ofstream outputFile(outputFileLocation);

  // obtain total icoing transitions for each src layer (used as denominator of
  // probability)
  std::vector<std::uint32_t> srcTotals;
  srcTotals.resize(m_detectorGeometry.size(), 0);

  for (const auto& [indexes, nTransitions] : m_layerPairs) {
    std::uint32_t layerIndex = getIndexByGbtsId(indexes.first);
    srcTotals[layerIndex] += nTransitions;
  }

  // find transitions that pass probability cut and add to temp container
  std::set<LayerIdPair> tempPairs;
  for (const auto& [indexes, nTransitions] : m_layerPairs) {
    const std::uint32_t srcIndex = getIndexByGbtsId(indexes.first);

    float probability{};

    if (srcTotals[srcIndex] == 0) {
      // avoid NAN error with 0/0
      probability = 0;
    } else {
      probability = static_cast<float>(nTransitions) / srcTotals[srcIndex];
    }

    const bool passCut = (m_cfg.probThreshold == -1)
                             ? (probability != 0)
                             : (probability >= m_cfg.probThreshold);

    if (passCut) {
      tempPairs.emplace(indexes.first, indexes.second);
    }
  }

  // if symmetrizing connection table, add transitions that mirror found ones
  if (m_cfg.doSymmetrization) {
    for (const auto& pair : tempPairs) {
      // find swapped ids
      const auto srcSwappedId = oppositeSideLayer(pair.first);
      const auto dstSwappedId = oppositeSideLayer(pair.second);

      if (!srcSwappedId && !dstSwappedId) {
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

  // finally, add transitions to output file (old or new format)
  if (m_cfg.useOldFormatting) {
    oldStyleFormatting(outputFileLocation, tempPairs);
  } else {
    for (const auto& pair : tempPairs) {
      // swap order as we want outward -> inward ordering
      outputFile << pair.second << " " << pair.first << "\n";
    }
  }
}

std::optional<std::int32_t> GbtsLayerConnectionTool::findGbtsIdByCoord(
    const float r, const float z) const {
  for (const auto& layer : m_detectorGeometry) {
    const float zMin = layer.minZ - m_cfg.zMinTol;
    const float zMax = layer.maxZ + m_cfg.zMaxTol;
    const float rMin = layer.minR - m_cfg.rMinTol;
    const float rMax = layer.maxR + m_cfg.rMaxTol;

    if (zMin <= z && z <= zMax) {
      if (rMin <= r && r <= rMax) {
        return layer.gbtsId;
      }
    }
  }

  return std::nullopt;
}

std::uint32_t GbtsLayerConnectionTool::getIndexByGbtsId(
    std::int32_t gbtsId) const {
  for (std::uint32_t idx = 0; idx < m_detectorGeometry.size(); idx++) {
    if (gbtsId == m_detectorGeometry[idx].gbtsId) {
      return idx;
    }
  }

  throw std::runtime_error("index not found for GBTS ID");
}

std::optional<std::int32_t> GbtsLayerConnectionTool::oppositeSideLayer(
    std::int32_t layerId) const {
  const std::uint32_t layerIndex = getIndexByGbtsId(layerId);

  const auto& layer = m_detectorGeometry[layerIndex];

  bool switchedMinZ{};
  bool switchedMaxZ{};

  bool sameMaxR{};
  bool sameMinR{};

  for (const auto& switchedLayer : m_detectorGeometry) {
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

void GbtsLayerConnectionTool::oldStyleFormatting(
    const std::string& outputFileLocations,
    const std::set<LayerIdPair>& tempTable) const {
  std::ofstream outputFile(outputFileLocations);

  outputFile << tempTable.size() << " " << 0.2 << "\n";
  for (const auto& pair : tempTable) {
    outputFile << 0 << " " << 1 << " " << pair.second << " " << pair.first
               << " " << 1 << " " << 1 << " " << 100 << "\n";

    outputFile << 100 << "\n";
  }
}

}  // namespace Acts::Experimental
