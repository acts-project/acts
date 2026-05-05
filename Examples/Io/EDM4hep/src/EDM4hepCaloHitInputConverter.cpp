// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepCaloHitInputConverter.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include <edm4hep/CaloHitContribution.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

namespace ActsExamples {

namespace {

constexpr std::uint64_t kUnmatched =
    std::numeric_limits<std::uint64_t>::max();

}  // namespace

EDM4hepCaloHitInputConverter::DetectorEncoder
EDM4hepCaloHitInputConverter::defaultDetectorEncoder() {
  // Codes match utils/detector_enums.py::CALO_DETECTOR_CODES.
  return [](std::string_view name, double z) -> std::uint8_t {
    if (name == "ECalBarrelCollection") {
      return 10;
    }
    if (name == "ECalEndcapCollection") {
      return z < 0.0 ? 9 : 11;
    }
    if (name == "HCalBarrelCollection") {
      return 13;
    }
    if (name == "HCalEndcapCollection") {
      return z < 0.0 ? 12 : 14;
    }
    return 255;
  };
}

EDM4hepCaloHitInputConverter::EDM4hepCaloHitInputConverter(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : PodioInputConverter("EDM4hepCaloHitInputConverter", cfg.inputFrame,
                          std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.inputCaloHitCollections.empty()) {
    throw std::invalid_argument("Missing calorimeter hit collections");
  }
  if (m_cfg.inputMCParticleMap.empty()) {
    throw std::invalid_argument("Missing MCParticle index map input");
  }
  if (m_cfg.outputCaloHits.empty()) {
    throw std::invalid_argument("Missing output calo hits container name");
  }
  if (!m_cfg.detectorEncoder) {
    throw std::invalid_argument("detectorEncoder must be set");
  }
  if (!m_cfg.isEcalCollection) {
    throw std::invalid_argument("isEcalCollection must be set");
  }
  if (m_cfg.lightSpeed <= 0.0) {
    throw std::invalid_argument("lightSpeed must be positive");
  }

  m_inputMCParticleMap.initialize(m_cfg.inputMCParticleMap);
  m_outputCaloHits.initialize(m_cfg.outputCaloHits);
}

ProcessCode EDM4hepCaloHitInputConverter::convert(
    const AlgorithmContext& ctx, const podio::Frame& frame) const {
  const EDM4hepMCParticleIndexMap& mcpMap = m_inputMCParticleMap(ctx);

  CaloHitContainer outCells;

  // (collection index, cellID) → index in `outCells`. Keying on the
  // collection too avoids cellID collisions across subdetectors.
  using CellKey = std::pair<std::size_t, std::uint64_t>;
  struct CellKeyHash {
    std::size_t operator()(const CellKey& k) const noexcept {
      return std::hash<std::uint64_t>{}(k.second) ^
             (std::hash<std::size_t>{}(k.first) << 1);
    }
  };
  std::unordered_map<CellKey, std::size_t, CellKeyHash> cellToIdx;

  for (std::size_t collIdx = 0; collIdx < m_cfg.inputCaloHitCollections.size();
       ++collIdx) {
    const std::string& collName = m_cfg.inputCaloHitCollections[collIdx];
    const auto* base = frame.get(collName);
    if (base == nullptr) {
      throw std::runtime_error(
          "EDM4hepCaloHitInputConverter: collection '" + collName +
          "' not present in frame");
    }
    const auto* hits =
        dynamic_cast<const edm4hep::SimCalorimeterHitCollection*>(base);
    if (hits == nullptr) {
      throw std::runtime_error(
          "EDM4hepCaloHitInputConverter: collection '" + collName +
          "' is not a SimCalorimeterHitCollection");
    }

    const bool isEcal = m_cfg.isEcalCollection(collName);
    const double timeMin = isEcal ? m_cfg.ecalTimeMin : m_cfg.hcalTimeMin;
    const double timeMax = isEcal ? m_cfg.ecalTimeMax : m_cfg.hcalTimeMax;

    for (const auto& hit : *hits) {
      const std::uint64_t cellId = hit.getCellID();
      const auto& pos = hit.getPosition();
      const CellKey key{collIdx, cellId};

      // Defer cell creation: if every contribution is dropped by the time
      // window we don't want to emit an empty cell.
      auto cellIt = cellToIdx.end();
      // Per-cell particle bucket; `byParticleIdx` maps the resolved
      // particleRow back into `cell.contributions` for O(1) updates while
      // preserving first-seen ordering.
      std::unordered_map<std::uint64_t, std::size_t> byParticleIdx;

      for (const auto& contrib : hit.getContributions()) {
        const auto& step = contrib.getStepPosition();
        const double r =
            std::sqrt(step.x * step.x + step.y * step.y + step.z * step.z);
        const double correctedTime =
            contrib.getTime() - (r / m_cfg.lightSpeed - m_cfg.tofOffset);
        if (correctedTime < timeMin || correctedTime > timeMax) {
          continue;
        }

        if (cellIt == cellToIdx.end()) {
          cellIt = cellToIdx.find(key);
          if (cellIt == cellToIdx.end()) {
            CaloHit cell;
            cell.cellId = cellId;
            cell.detector = m_cfg.detectorEncoder(collName, pos.z);
            cell.position = Acts::Vector3{pos.x, pos.y, pos.z};
            cell.totalEnergy = 0.0F;
            outCells.push_back(std::move(cell));
            cellIt = cellToIdx.emplace(key, outCells.size() - 1).first;
          } else {
            // Pre-populate from the stored contributions so subsequent
            // contributions to the same particle aggregate correctly. This
            // path is reachable only if the same (collection, cellID) shows
            // up more than once in the input collection, which shouldn't
            // happen with well-formed EDM4hep but we keep the code honest.
            const CaloHit& existing = outCells[cellIt->second];
            byParticleIdx.reserve(existing.contributions.size());
            for (std::size_t i = 0; i < existing.contributions.size(); ++i) {
              byParticleIdx.emplace(existing.contributions[i].particleRow, i);
            }
          }
        }
        CaloHit& cell = outCells[cellIt->second];

        const auto particle = contrib.getParticle();
        std::uint64_t particleRow = kUnmatched;
        if (particle.isAvailable()) {
          auto mIt = mcpMap.find(particle.getObjectID().index);
          if (mIt != mcpMap.end()) {
            particleRow = static_cast<std::uint64_t>(mIt->second);
          }
        }

        const float e = static_cast<float>(contrib.getEnergy());
        const float t = contrib.getTime();

        auto pIt = byParticleIdx.find(particleRow);
        if (pIt == byParticleIdx.end()) {
          // Store sums so we can compute the energy-weighted time at the
          // end. We re-use CaloHitContribution::time as a running Σ(e·t)
          // here and divide by Σe in the final pass.
          CaloHitContribution acc;
          acc.particleRow = particleRow;
          acc.energy = e;
          acc.time = e * t;
          cell.contributions.push_back(acc);
          byParticleIdx.emplace(particleRow, cell.contributions.size() - 1);
        } else {
          CaloHitContribution& acc = cell.contributions[pIt->second];
          acc.energy += e;
          acc.time += e * t;  // running Σ(e·t); finalised below
        }
        cell.totalEnergy += e;
      }
    }
  }

  // Finalise: convert the stored Σ(e·t) into the energy-weighted average.
  for (CaloHit& cell : outCells) {
    for (CaloHitContribution& acc : cell.contributions) {
      acc.time = acc.energy != 0.0F ? acc.time / acc.energy : 0.0F;
    }
  }

  ACTS_DEBUG("Produced " << outCells.size() << " calorimeter cells");
  m_outputCaloHits(ctx, std::move(outCells));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
