// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepCaloHitInputConverter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <edm4hep/CaloHitContribution.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

namespace ActsExamples {

namespace {

constexpr std::uint64_t kUnmatched = std::numeric_limits<std::uint64_t>::max();

}  // namespace

CaloCollectionDetectorCodes CaloCollectionDetectorCodes::barrel(
    std::uint8_t code) {
  CaloCollectionDetectorCodes out{};
  out.isBarrel = true;
  out.barrelCode = code;
  return out;
}

CaloCollectionDetectorCodes CaloCollectionDetectorCodes::endcap(
    std::uint8_t negZ, std::uint8_t posZ) {
  CaloCollectionDetectorCodes out{};
  out.isBarrel = false;
  out.endcapNegCode = negZ;
  out.endcapPosCode = posZ;
  return out;
}

std::vector<CaloCollectionDetectorCodes>
EDM4hepCaloHitInputConverter::buildDetectorCodesByCollectionIndex(
    const std::vector<std::string>& inputCaloHitCollections,
    const std::unordered_map<std::string, CaloCollectionDetectorCodes>&
        caloDetectorCodesByCollectionName) {
  static const CaloCollectionDetectorCodes kUnknown =
      CaloCollectionDetectorCodes::barrel(255);
  std::vector<CaloCollectionDetectorCodes> out;
  out.reserve(inputCaloHitCollections.size());
  for (const std::string& collName : inputCaloHitCollections) {
    const auto it = caloDetectorCodesByCollectionName.find(collName);
    out.push_back(it != caloDetectorCodesByCollectionName.end() ? it->second
                                                                : kUnknown);
  }
  return out;
}

std::function<bool(std::string_view)>
EDM4hepCaloHitInputConverter::defaultIsEcalCollection() {
  return [](std::string_view name) {
    if (name.starts_with("ECal")) {
      return true;
    }
    if (name.starts_with("HCal")) {
      return false;
    }
    throw std::invalid_argument(
        std::string("EDM4hepCaloHitInputConverter: cannot infer "
                    "ECal/HCal from collection name '") +
        std::string(name) + "'");
  };
}

EDM4hepCaloHitInputConverter::EDM4hepCaloHitInputConverter(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : PodioInputConverter("EDM4hepCaloHitInputConverter", cfg.inputFrame,
                          std::move(logger)),
      m_cfg(cfg),
      m_detectorCodesByCollectionIndex(buildDetectorCodesByCollectionIndex(
          m_cfg.inputCaloHitCollections,
          m_cfg.caloDetectorCodesByCollectionName)) {
  if (m_cfg.inputCaloHitCollections.empty()) {
    throw std::invalid_argument("Missing calorimeter hit collections");
  }
  if (m_cfg.inputMCParticleMap.empty()) {
    throw std::invalid_argument("Missing MCParticle index map input");
  }
  if (m_cfg.outputCaloHits.empty()) {
    throw std::invalid_argument("Missing output calo hits container name");
  }
  if (!m_cfg.isEcalCollection) {
    throw std::invalid_argument("isEcalCollection must be set");
  }

  m_inputMCParticleMap.initialize(m_cfg.inputMCParticleMap);
  m_outputCaloHits.initialize(m_cfg.outputCaloHits);

  // Section-level averaging timers shared across the whole job. One sample
  // per event per section; logs aggregate stats at algorithm teardown.
  m_timerResolve.emplace("Resolving collections", this->logger(),
                         Acts::Logging::INFO);
  m_timerAggregate.emplace("Aggregating contributions", this->logger(),
                           Acts::Logging::INFO);
  m_timerFinalise.emplace("Finalising contribution times", this->logger(),
                          Acts::Logging::INFO);
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

  // Per-collection accounting so log output can pinpoint which subdetector
  // is losing cells (timing-window cuts, missing input, encoder fall-through,
  // etc.).
  struct CollStats {
    std::size_t nHits = 0;
    std::size_t nContribsTotal = 0;
    std::size_t nContribsKept = 0;
    std::size_t nCellsEmitted = 0;
  };
  std::vector<CollStats> stats(m_cfg.inputCaloHitCollections.size());

  // Resolve & validate all input collections up front so we can size the
  // output container and the cell-lookup map to the worst case (one cell
  // per input hit). Reserving keeps the hot loop free of vector
  // reallocations and unordered_map rehashes for million-hit events.
  std::vector<const edm4hep::SimCalorimeterHitCollection*> hitColls(
      m_cfg.inputCaloHitCollections.size(), nullptr);
  std::size_t totalHits = 0;
  {
    auto resolveSample = m_timerResolve->sample();
    for (std::size_t collIdx = 0;
         collIdx < m_cfg.inputCaloHitCollections.size(); ++collIdx) {
      const std::string& collName = m_cfg.inputCaloHitCollections[collIdx];
      const auto* base = frame.get(collName);
      if (base == nullptr) {
        throw std::runtime_error("EDM4hepCaloHitInputConverter: collection '" +
                                 collName + "' not present in frame");
      }
      const auto* hits =
          dynamic_cast<const edm4hep::SimCalorimeterHitCollection*>(base);
      if (hits == nullptr) {
        throw std::runtime_error("EDM4hepCaloHitInputConverter: collection '" +
                                 collName +
                                 "' is not a SimCalorimeterHitCollection");
      }
      hitColls[collIdx] = hits;
      totalHits += hits->size();
      stats[collIdx].nHits = hits->size();
    }
    outCells.reserve(totalHits);
    cellToIdx.reserve(totalHits);
  }

  // Per-cell particle bucket; `byParticleIdx` maps the resolved
  // `particleRow` back into `cell.contributions` for O(1) updates while
  // preserving first-seen ordering. Lifted out of the hit loop so its
  // bucket array is reused across millions of hits — `clear()` keeps the
  // existing capacity, allocating only on first growth.
  std::unordered_map<std::uint64_t, std::size_t> byParticleIdx;

  {
    auto aggregateSample = m_timerAggregate->sample();

    for (std::size_t collIdx = 0; collIdx < hitColls.size(); ++collIdx) {
      const std::string& collName = m_cfg.inputCaloHitCollections[collIdx];
      const auto* hits = hitColls[collIdx];

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
        byParticleIdx.clear();

        // The TOF radius is taken from the cell centre, matching pyedm4hep's
        // `EDM4hepEventBatch._load_calo_contributions`, which merges the
        // cell (x, y, z) onto every contribution row before computing
        // `r = sqrt(x²+y²+z²)`. The per-contribution `stepPosition` is left
        // unused on purpose; many sims (e.g. the default Geant4 calo SD
        // used by ODD) don't even fill it. EDM4hep positions are in mm,
        // which is the native ACTS length unit, so no conversion is needed.
        const double r =
            std::sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
        const double tofShift =
            r / Acts::PhysicalConstants::c - m_cfg.tofOffset;

        for (const auto& contrib : hit.getContributions()) {
          stats[collIdx].nContribsTotal += 1;
          // EDM4hep contribution times are in ns; convert to ACTS time
          // units so the comparison with `timeMin`/`timeMax` and the
          // tofShift (both in ACTS units) stays consistent.
          const double correctedTime =
              contrib.getTime() * Acts::UnitConstants::ns - tofShift;
          if (correctedTime < timeMin || correctedTime > timeMax) {
            continue;
          }
          stats[collIdx].nContribsKept += 1;

          if (cellIt == cellToIdx.end()) {
            cellIt = cellToIdx.find(key);
            if (cellIt == cellToIdx.end()) {
              CaloHit cell;
              cell.cellId = cellId;
              const CaloCollectionDetectorCodes& detCodes =
                  m_detectorCodesByCollectionIndex[collIdx];
              cell.detector = detCodes.detectorCode(pos.z);
              cell.position = Acts::Vector3{pos.x, pos.y, pos.z};
              cell.totalEnergy = 0.0F;
              outCells.push_back(std::move(cell));
              stats[collIdx].nCellsEmitted += 1;
              cellIt = cellToIdx.emplace(key, outCells.size() - 1).first;
            } else {
              // Pre-populate from the stored contributions so subsequent
              // contributions to the same particle aggregate correctly.
              // This path is reachable only if the same (collection,
              // cellID) shows up more than once in the input collection,
              // which shouldn't happen with well-formed EDM4hep but we keep
              // the code honest.
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

          const float e = contrib.getEnergy();
          const float t = contrib.getTime();

          auto pIt = byParticleIdx.find(particleRow);
          if (pIt == byParticleIdx.end()) {
            // Store sums so we can compute the energy-weighted time at the
            // end. We reuse CaloHitContribution::time as a running Σ(e·t)
            // here and divide by Σe in the final pass.
            CaloHitContribution acc{};
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
  }

  // Finalise: convert the stored Σ(e·t) into the energy-weighted average.
  {
    auto finaliseSample = m_timerFinalise->sample();
    for (CaloHit& cell : outCells) {
      for (CaloHitContribution& acc : cell.contributions) {
        acc.time = acc.energy != 0.0F ? acc.time / acc.energy : 0.0F;
      }
    }
  }

  for (std::size_t i = 0; i < m_cfg.inputCaloHitCollections.size(); ++i) {
    ACTS_DEBUG("calo collection '"
               << m_cfg.inputCaloHitCollections[i] << "': hits="
               << stats[i].nHits << " contribs=" << stats[i].nContribsTotal
               << " kept=" << stats[i].nContribsKept
               << " cells_emitted=" << stats[i].nCellsEmitted);
  }

  ACTS_DEBUG("Produced " << outCells.size() << " calorimeter cells");
  m_outputCaloHits(ctx, std::move(outCells));
  return ProcessCode::SUCCESS;
}

ProcessCode EDM4hepCaloHitInputConverter::finalize() {
  // Resetting triggers output
  m_timerResolve.reset();
  m_timerAggregate.reset();
  m_timerFinalise.reset();
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
