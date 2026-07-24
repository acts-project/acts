// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Arrow/ArrowCaloHitOutputConverter.hpp"

#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <arrow/api.h>

namespace ActsExamples {

namespace {

void check(const arrow::Status& s, const char* what) {
  if (!s.ok()) {
    throw std::runtime_error(std::string(what) + ": " + s.ToString());
  }
}

}  // namespace

ArrowCaloHitOutputConverter::CellThresholdFn
ArrowCaloHitOutputConverter::defaultCellThreshold(double ecalEnergyThreshold,
                                                  double hcalEnergyThreshold) {
  return
      [ecalEnergyThreshold, hcalEnergyThreshold](std::uint8_t det) -> double {
        if (det >= 9 && det <= 11) {
          return ecalEnergyThreshold;
        }
        if (det >= 12 && det <= 14) {
          return hcalEnergyThreshold;
        }
        return 0.0;
      };
}

ArrowCaloHitOutputConverter::ArrowCaloHitOutputConverter(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : ArrowOutputConverter("ArrowCaloHitOutputConverter", std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.inputCaloHits.empty()) {
    throw std::invalid_argument("Missing calo hits input collection");
  }
  if (m_cfg.outputTable.empty()) {
    throw std::invalid_argument("Missing output table name");
  }

  m_cellThreshold = m_cfg.cellThreshold
                        ? m_cfg.cellThreshold
                        : defaultCellThreshold(m_cfg.ecalEnergyThreshold,
                                               m_cfg.hcalEnergyThreshold);

  m_inputCaloHits.initialize(m_cfg.inputCaloHits);
  m_outputTable.initialize(m_cfg.outputTable);
}

std::vector<std::string> ArrowCaloHitOutputConverter::collections() const {
  return {m_cfg.outputTable};
}

ProcessCode ArrowCaloHitOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  const CaloHitContainer& cells = m_inputCaloHits(ctx);

  auto* pool = arrow::default_memory_pool();

  arrow::ListBuilder detList(pool, std::make_shared<arrow::UInt8Builder>(pool));
  arrow::ListBuilder energyList(pool,
                                std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder xList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder yList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder zList(pool, std::make_shared<arrow::FloatBuilder>(pool));

  // Two-level nested builders: the outer list has one element per event row
  // (we open it once below), and each inner list is one cell's contribution
  // values.
  auto pidInner = std::make_shared<arrow::ListBuilder>(
      pool, std::make_shared<arrow::UInt64Builder>(pool));
  auto eInner = std::make_shared<arrow::ListBuilder>(
      pool, std::make_shared<arrow::FloatBuilder>(pool));
  auto tInner = std::make_shared<arrow::ListBuilder>(
      pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pidList(pool, pidInner);
  arrow::ListBuilder eList(pool, eInner);
  arrow::ListBuilder tList(pool, tInner);

  check(detList.Append(), "open detector list");
  check(energyList.Append(), "open total_energy list");
  check(xList.Append(), "open x list");
  check(yList.Append(), "open y list");
  check(zList.Append(), "open z list");
  check(pidList.Append(), "open contrib_particle_ids list");
  check(eList.Append(), "open contrib_energies list");
  check(tList.Append(), "open contrib_times list");

  auto* detV = static_cast<arrow::UInt8Builder*>(detList.value_builder());
  auto* energyV = static_cast<arrow::FloatBuilder*>(energyList.value_builder());
  auto* xV = static_cast<arrow::FloatBuilder*>(xList.value_builder());
  auto* yV = static_cast<arrow::FloatBuilder*>(yList.value_builder());
  auto* zV = static_cast<arrow::FloatBuilder*>(zList.value_builder());

  auto* pidVL = static_cast<arrow::ListBuilder*>(pidList.value_builder());
  auto* eVL = static_cast<arrow::ListBuilder*>(eList.value_builder());
  auto* tVL = static_cast<arrow::ListBuilder*>(tList.value_builder());
  auto* pidVV = static_cast<arrow::UInt64Builder*>(pidVL->value_builder());
  auto* eVV = static_cast<arrow::FloatBuilder*>(eVL->value_builder());
  auto* tVV = static_cast<arrow::FloatBuilder*>(tVL->value_builder());

  // Conservative reservation; at most every cell survives.
  const std::size_t nCells = cells.size();
  check(detV->Reserve(nCells), "reserve detector");
  check(energyV->Reserve(nCells), "reserve total_energy");
  check(xV->Reserve(nCells), "reserve x");
  check(yV->Reserve(nCells), "reserve y");
  check(zV->Reserve(nCells), "reserve z");

  for (const CaloHit& cell : cells) {
    const double threshold = m_cellThreshold(cell.detector);
    if (cell.totalEnergy < threshold) {
      continue;
    }

    detV->UnsafeAppend(cell.detector);
    energyV->UnsafeAppend(cell.totalEnergy);
    xV->UnsafeAppend(static_cast<float>(cell.position.x()));
    yV->UnsafeAppend(static_cast<float>(cell.position.y()));
    zV->UnsafeAppend(static_cast<float>(cell.position.z()));

    check(pidVL->Append(), "open per-cell particle list");
    check(eVL->Append(), "open per-cell energy list");
    check(tVL->Append(), "open per-cell time list");

    const std::size_t nContribs = cell.contributions.size();
    check(pidVV->Reserve(nContribs), "reserve cell particle ids");
    check(eVV->Reserve(nContribs), "reserve cell contrib energies");
    check(tVV->Reserve(nContribs), "reserve cell contrib times");
    for (const CaloHitContribution& acc : cell.contributions) {
      pidVV->UnsafeAppend(acc.particleRow);
      eVV->UnsafeAppend(acc.energy);
      tVV->UnsafeAppend(acc.time);
    }
  }

  auto finish = [](arrow::ListBuilder& b) {
    std::shared_ptr<arrow::Array> out;
    check(b.Finish(&out), "finish list");
    return out;
  };

  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      finish(detList), finish(energyList), finish(xList), finish(yList),
      finish(zList),   finish(pidList),    finish(eList), finish(tList),
  };

  auto table =
      arrow::Table::Make(ActsPlugins::ArrowUtil::caloHitSchema(), arrays);
  m_outputTable(ctx, ActsPlugins::ArrowUtil::ArrowTable{std::move(table)});

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
