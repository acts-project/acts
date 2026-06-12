// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Arrow/ArrowSimHitOutputConverter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>
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

ArrowSimHitOutputConverter::ArrowSimHitOutputConverter(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : ArrowOutputConverter("ArrowSimHitOutputConverter", std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing sim hits input collection");
  }
  if (m_cfg.outputTable.empty()) {
    throw std::invalid_argument("Missing output table name");
  }
  if (!m_cfg.detectorResolver) {
    throw std::invalid_argument("detectorResolver must be set");
  }

  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputParticles.maybeInitialize(m_cfg.inputParticles);
  m_outputTable.initialize(m_cfg.outputTable);
}

std::function<std::uint8_t(Acts::GeometryIdentifier)>
ArrowSimHitOutputConverter::makeVolumeIdDetectorResolver(
    const std::unordered_map<std::uint32_t, std::uint8_t>& volumeToDetector,
    std::uint8_t defaultValue) {
  constexpr std::size_t kNumVolumes =
      static_cast<std::size_t>(Acts::GeometryIdentifier::getMaxVolume()) + 1;
  std::array<std::uint8_t, kNumVolumes> detectorArray{};
  detectorArray.fill(defaultValue);
  for (const auto& [volume, detector] : volumeToDetector) {
    if (volume >= kNumVolumes) {
      throw std::invalid_argument(
          "makeVolumeIdDetectorResolver: volume id " + std::to_string(volume) +
          " exceeds maximum " +
          std::to_string(Acts::GeometryIdentifier::getMaxVolume()));
    }
    detectorArray[volume] = detector;
  }
  return [detectorArray](Acts::GeometryIdentifier gid) -> std::uint8_t {
    const auto volume = static_cast<std::size_t>(gid.volume());
    return detectorArray[volume];
  };
}

std::vector<std::string> ArrowSimHitOutputConverter::collections() const {
  return {m_cfg.outputTable};
}

ProcessCode ArrowSimHitOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  const SimHitContainer& simHits = m_inputSimHits(ctx);
  const SimParticleContainer* particles =
      m_inputParticles.isInitialized() ? &m_inputParticles(ctx) : nullptr;

  auto* pool = arrow::default_memory_pool();

  // Truth table: one entry per sim-hit, ALL sim-hits in container order. The
  // entry's position is the sim-hit id that the measurement table references
  // via simhit_ids.
  arrow::ListBuilder txList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder tyList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder tzList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder ttList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder tpxList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder tpyList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder tpzList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder teList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder deList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pidList(pool,
                             std::make_shared<arrow::UInt64Builder>(pool));
  arrow::ListBuilder hidxList(pool,
                              std::make_shared<arrow::UInt16Builder>(pool));
  arrow::ListBuilder detList(pool, std::make_shared<arrow::UInt8Builder>(pool));
  arrow::ListBuilder volList(pool, std::make_shared<arrow::UInt8Builder>(pool));
  arrow::ListBuilder layList(pool,
                             std::make_shared<arrow::UInt16Builder>(pool));
  arrow::ListBuilder surfList(pool,
                              std::make_shared<arrow::UInt32Builder>(pool));

  std::array<arrow::ListBuilder*, 15> all = {
      &txList,  &tyList,  &tzList, &ttList,   &tpxList,
      &tpyList, &tpzList, &teList, &deList,   &pidList,
      &hidxList, &detList, &volList, &layList, &surfList};
  for (auto* b : all) {
    check(b->Append(), "open per-event list");
  }

  auto* txV = static_cast<arrow::FloatBuilder*>(txList.value_builder());
  auto* tyV = static_cast<arrow::FloatBuilder*>(tyList.value_builder());
  auto* tzV = static_cast<arrow::FloatBuilder*>(tzList.value_builder());
  auto* ttV = static_cast<arrow::FloatBuilder*>(ttList.value_builder());
  auto* tpxV = static_cast<arrow::FloatBuilder*>(tpxList.value_builder());
  auto* tpyV = static_cast<arrow::FloatBuilder*>(tpyList.value_builder());
  auto* tpzV = static_cast<arrow::FloatBuilder*>(tpzList.value_builder());
  auto* teV = static_cast<arrow::FloatBuilder*>(teList.value_builder());
  auto* deV = static_cast<arrow::FloatBuilder*>(deList.value_builder());
  auto* pidV = static_cast<arrow::UInt64Builder*>(pidList.value_builder());
  auto* hidxV = static_cast<arrow::UInt16Builder*>(hidxList.value_builder());
  auto* detV = static_cast<arrow::UInt8Builder*>(detList.value_builder());
  auto* volV = static_cast<arrow::UInt8Builder*>(volList.value_builder());
  auto* layV = static_cast<arrow::UInt16Builder*>(layList.value_builder());
  auto* surfV = static_cast<arrow::UInt32Builder*>(surfList.value_builder());

  const std::size_t nHits = simHits.size();
  check(txV->Reserve(nHits), "reserve true_x");
  check(tyV->Reserve(nHits), "reserve true_y");
  check(tzV->Reserve(nHits), "reserve true_z");
  check(ttV->Reserve(nHits), "reserve true_time");
  check(tpxV->Reserve(nHits), "reserve tpx");
  check(tpyV->Reserve(nHits), "reserve tpy");
  check(tpzV->Reserve(nHits), "reserve tpz");
  check(teV->Reserve(nHits), "reserve tE");
  check(deV->Reserve(nHits), "reserve dE");
  check(pidV->Reserve(nHits), "reserve particle_id");
  check(hidxV->Reserve(nHits), "reserve hit_index");
  check(detV->Reserve(nHits), "reserve detector");
  check(volV->Reserve(nHits), "reserve volume_id");
  check(layV->Reserve(nHits), "reserve layer_id");
  check(surfV->Reserve(nHits), "reserve surface_id");

  constexpr std::uint64_t kUnmatched =
      std::numeric_limits<std::uint64_t>::max();
  constexpr std::uint16_t kNoIndex = std::numeric_limits<std::uint16_t>::max();

  for (const auto& hit : simHits) {
    const auto& pos = hit.fourPosition();
    txV->UnsafeAppend(static_cast<float>(pos.x() / Acts::UnitConstants::mm));
    tyV->UnsafeAppend(static_cast<float>(pos.y() / Acts::UnitConstants::mm));
    tzV->UnsafeAppend(static_cast<float>(pos.z() / Acts::UnitConstants::mm));
    ttV->UnsafeAppend(static_cast<float>(pos.w() / Acts::UnitConstants::mm));

    const auto& mom = hit.momentum4Before();
    tpxV->UnsafeAppend(static_cast<float>(mom.x() / Acts::UnitConstants::GeV));
    tpyV->UnsafeAppend(static_cast<float>(mom.y() / Acts::UnitConstants::GeV));
    tpzV->UnsafeAppend(static_cast<float>(mom.z() / Acts::UnitConstants::GeV));
    teV->UnsafeAppend(static_cast<float>(mom.w() / Acts::UnitConstants::GeV));
    deV->UnsafeAppend(
        static_cast<float>(hit.depositedEnergy() / Acts::UnitConstants::GeV));

    std::uint64_t pid = kUnmatched;
    if (particles != nullptr) {
      auto pIt = particles->find(hit.particleId());
      if (pIt != particles->end()) {
        pid =
            static_cast<std::uint64_t>(std::distance(particles->begin(), pIt));
      }
    }
    pidV->UnsafeAppend(pid);

    // Hit index along the particle trajectory (ordering for consumers); -1
    // (unset) maps to the uint16 sentinel.
    const std::int32_t idx = hit.index();
    hidxV->UnsafeAppend(
        (idx < 0 || idx > static_cast<std::int32_t>(kNoIndex))
            ? kNoIndex
            : static_cast<std::uint16_t>(idx));

    const auto gid = hit.geometryId();
    detV->UnsafeAppend(m_cfg.detectorResolver(gid));
    volV->UnsafeAppend(static_cast<std::uint8_t>(gid.volume()));
    layV->UnsafeAppend(static_cast<std::uint16_t>(gid.layer()));
    surfV->UnsafeAppend(static_cast<std::uint32_t>(gid.sensitive()));
  }

  auto finish = [](arrow::ListBuilder& b) {
    std::shared_ptr<arrow::Array> out;
    check(b.Finish(&out), "finish list");
    return out;
  };

  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      finish(txList),  finish(tyList),  finish(tzList),  finish(ttList),
      finish(tpxList), finish(tpyList), finish(tpzList), finish(teList),
      finish(deList),  finish(pidList), finish(hidxList), finish(detList),
      finish(volList), finish(layList), finish(surfList),
  };

  auto table =
      arrow::Table::Make(ActsPlugins::ArrowUtil::simHitSchema(), arrays);
  m_outputTable(ctx, ActsPlugins::ArrowUtil::ArrowTable{std::move(table)});

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
