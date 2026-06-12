// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Arrow/ArrowSimHitOutputConverter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <array>
#include <cmath>
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

  // v2 emits one row per measurement, so the measurement container, the
  // sim-hit -> measurement map (inverted internally), and the tracking
  // geometry (for the reco position) are all required.
  if (m_cfg.inputMeasurements.empty() ||
      m_cfg.inputSimHitMeasurementsMap.empty() ||
      m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument(
        "ArrowSimHitOutputConverter: inputMeasurements, "
        "inputSimHitMeasurementsMap, and trackingGeometry are required (one "
        "row per measurement)");
  }

  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputParticles.maybeInitialize(m_cfg.inputParticles);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSimHitMeasurementsMap.initialize(m_cfg.inputSimHitMeasurementsMap);
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
  const MeasurementContainer& measurements = m_inputMeasurements(ctx);
  const SimHitMeasurementsMap& simHitMeasMap = m_inputSimHitMeasurementsMap(ctx);

  // Invert the sim-hit -> measurement map into measurement -> sim-hits, and
  // record which sim-hits contribute to a measurement (the rest are dropped
  // from the v2 table and counted for the diagnostic below).
  std::unordered_multimap<Index, SimHitIndex> measToSimHits;
  measToSimHits.reserve(simHitMeasMap.size());
  std::vector<bool> referenced(simHits.size(), false);
  for (const auto& [simHitIdx, measIdx] : simHitMeasMap) {
    measToSimHits.emplace(measIdx, simHitIdx);
    if (simHitIdx < referenced.size()) {
      referenced[simHitIdx] = true;
    }
  }

  // Diagnostic: sim-hits belonging to no measurement are dropped from v2. This
  // is expected to be a tiny fraction (<0.01%); warn always and fail if it
  // exceeds the configured threshold, so a regression can't silently grow.
  const std::size_t totalSimHits = simHits.size();
  std::size_t droppedSimHits = 0;
  for (bool seen : referenced) {
    if (!seen) {
      ++droppedSimHits;
    }
  }
  const double droppedFraction =
      totalSimHits > 0
          ? static_cast<double>(droppedSimHits) / static_cast<double>(totalSimHits)
          : 0.0;
  if (droppedSimHits > 0) {
    ACTS_WARNING("ArrowSimHitOutputConverter: "
                 << droppedSimHits << " / " << totalSimHits << " sim-hits ("
                 << droppedFraction * 100.0
                 << "%) belong to no measurement and were dropped from the "
                    "measurement table");
  }
  if (droppedFraction > m_cfg.maxUnmatchedSimHitFraction) {
    throw std::runtime_error(
        "ArrowSimHitOutputConverter: dropped sim-hit fraction " +
        std::to_string(droppedFraction) + " exceeds threshold " +
        std::to_string(m_cfg.maxUnmatchedSimHitFraction));
  }

  auto* pool = arrow::default_memory_pool();

  // Flat columns: one value per measurement.
  arrow::ListBuilder xList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder yList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder zList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder detList(pool, std::make_shared<arrow::UInt8Builder>(pool));
  arrow::ListBuilder volList(pool, std::make_shared<arrow::UInt8Builder>(pool));
  arrow::ListBuilder layList(pool,
                             std::make_shared<arrow::UInt16Builder>(pool));
  arrow::ListBuilder surfList(pool,
                              std::make_shared<arrow::UInt32Builder>(pool));

  // Nested truth columns: outer = per measurement, inner = per contributing
  // sim-hit (mirrors the calo contrib_* columns).
  auto pidInner = std::make_shared<arrow::ListBuilder>(
      pool, std::make_shared<arrow::UInt64Builder>(pool));
  auto txInner = std::make_shared<arrow::ListBuilder>(
      pool, std::make_shared<arrow::FloatBuilder>(pool));
  auto tyInner = std::make_shared<arrow::ListBuilder>(
      pool, std::make_shared<arrow::FloatBuilder>(pool));
  auto tzInner = std::make_shared<arrow::ListBuilder>(
      pool, std::make_shared<arrow::FloatBuilder>(pool));
  auto timeInner = std::make_shared<arrow::ListBuilder>(
      pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pidList(pool, pidInner);
  arrow::ListBuilder txList(pool, txInner);
  arrow::ListBuilder tyList(pool, tyInner);
  arrow::ListBuilder tzList(pool, tzInner);
  arrow::ListBuilder timeList(pool, timeInner);

  // Open the per-event (outer) list on every column once.
  check(xList.Append(), "open x list");
  check(yList.Append(), "open y list");
  check(zList.Append(), "open z list");
  check(detList.Append(), "open detector list");
  check(volList.Append(), "open volume_id list");
  check(layList.Append(), "open layer_id list");
  check(surfList.Append(), "open surface_id list");
  check(pidList.Append(), "open particle_id list");
  check(txList.Append(), "open true_x list");
  check(tyList.Append(), "open true_y list");
  check(tzList.Append(), "open true_z list");
  check(timeList.Append(), "open time list");

  auto* xV = static_cast<arrow::FloatBuilder*>(xList.value_builder());
  auto* yV = static_cast<arrow::FloatBuilder*>(yList.value_builder());
  auto* zV = static_cast<arrow::FloatBuilder*>(zList.value_builder());
  auto* detV = static_cast<arrow::UInt8Builder*>(detList.value_builder());
  auto* volV = static_cast<arrow::UInt8Builder*>(volList.value_builder());
  auto* layV = static_cast<arrow::UInt16Builder*>(layList.value_builder());
  auto* surfV = static_cast<arrow::UInt32Builder*>(surfList.value_builder());

  auto* pidVL = static_cast<arrow::ListBuilder*>(pidList.value_builder());
  auto* txVL = static_cast<arrow::ListBuilder*>(txList.value_builder());
  auto* tyVL = static_cast<arrow::ListBuilder*>(tyList.value_builder());
  auto* tzVL = static_cast<arrow::ListBuilder*>(tzList.value_builder());
  auto* timeVL = static_cast<arrow::ListBuilder*>(timeList.value_builder());
  auto* pidVV = static_cast<arrow::UInt64Builder*>(pidVL->value_builder());
  auto* txVV = static_cast<arrow::FloatBuilder*>(txVL->value_builder());
  auto* tyVV = static_cast<arrow::FloatBuilder*>(tyVL->value_builder());
  auto* tzVV = static_cast<arrow::FloatBuilder*>(tzVL->value_builder());
  auto* timeVV = static_cast<arrow::FloatBuilder*>(timeVL->value_builder());

  const std::size_t nMeas = measurements.size();
  check(xV->Reserve(nMeas), "reserve x");
  check(yV->Reserve(nMeas), "reserve y");
  check(zV->Reserve(nMeas), "reserve z");
  check(detV->Reserve(nMeas), "reserve detector");
  check(volV->Reserve(nMeas), "reserve volume_id");
  check(layV->Reserve(nMeas), "reserve layer_id");
  check(surfV->Reserve(nMeas), "reserve surface_id");

  constexpr std::uint64_t kUnmatched =
      std::numeric_limits<std::uint64_t>::max();
  constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();

  for (Index measIdx = 0; measIdx < nMeas; ++measIdx) {
    const auto meas = measurements.getMeasurement(measIdx);
    const auto gid = meas.geometryId();

    // Reco (digitized) global position: project the measurement's bound
    // parameters through its surface. NaN for non-regular surfaces.
    float gx = kNaN;
    float gy = kNaN;
    float gz = kNaN;
    const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(gid);
    if (surface == nullptr) {
      throw std::runtime_error(
          "ArrowSimHitOutputConverter: surface not found for geometry id " +
          std::to_string(gid.value()));
    }
    if (const auto* regular =
            dynamic_cast<const Acts::RegularSurface*>(surface)) {
      // Expand the (possibly 1D) measurement to a full bound vector so local
      // positions are read from a stable layout regardless of the measured
      // subspace.
      const auto full = meas.fullParameters();
      Acts::Vector2 loc{full[Acts::eBoundLoc0], full[Acts::eBoundLoc1]};
      const Acts::Vector3 global = regular->localToGlobal(ctx.geoContext, loc);
      gx = static_cast<float>(global.x() / Acts::UnitConstants::mm);
      gy = static_cast<float>(global.y() / Acts::UnitConstants::mm);
      gz = static_cast<float>(global.z() / Acts::UnitConstants::mm);
    }
    xV->UnsafeAppend(gx);
    yV->UnsafeAppend(gy);
    zV->UnsafeAppend(gz);

    volV->UnsafeAppend(static_cast<std::uint8_t>(gid.volume()));
    layV->UnsafeAppend(static_cast<std::uint16_t>(gid.layer()));
    surfV->UnsafeAppend(static_cast<std::uint32_t>(gid.sensitive()));
    // The default `detectorResolver` reads the geometry id's `extra` byte,
    // which geometry construction is expected to stamp with a per-surface
    // subsystem id; users can swap in a custom resolver.
    detV->UnsafeAppend(m_cfg.detectorResolver(gid));

    // Open one inner list per measurement on each nested truth column.
    check(pidVL->Append(), "open per-measurement particle_id list");
    check(txVL->Append(), "open per-measurement true_x list");
    check(tyVL->Append(), "open per-measurement true_y list");
    check(tzVL->Append(), "open per-measurement true_z list");
    check(timeVL->Append(), "open per-measurement time list");

    // Reserve the nested inner value builders for this measurement's
    // contributors, then append. (Unlike the flat builders we cannot reserve
    // these up front because the per-measurement multiplicity is only known
    // here; without the reserve, UnsafeAppend would write past the buffer.)
    auto range = measToSimHits.equal_range(measIdx);
    const auto nContribs =
        static_cast<std::int64_t>(std::distance(range.first, range.second));
    check(pidVV->Reserve(nContribs), "reserve particle_id contribs");
    check(txVV->Reserve(nContribs), "reserve true_x contribs");
    check(tyVV->Reserve(nContribs), "reserve true_y contribs");
    check(tzVV->Reserve(nContribs), "reserve true_z contribs");
    check(timeVV->Reserve(nContribs), "reserve time contribs");
    for (auto it = range.first; it != range.second; ++it) {
      const SimHitIndex simHitIdx = it->second;
      const auto& hit = *simHits.nth(simHitIdx);
      const auto& pos = hit.fourPosition();
      txVV->UnsafeAppend(static_cast<float>(pos.x() / Acts::UnitConstants::mm));
      tyVV->UnsafeAppend(static_cast<float>(pos.y() / Acts::UnitConstants::mm));
      tzVV->UnsafeAppend(static_cast<float>(pos.z() / Acts::UnitConstants::mm));
      timeVV->UnsafeAppend(
          static_cast<float>(pos.w() / Acts::UnitConstants::mm));

      std::uint64_t pid = kUnmatched;
      if (particles != nullptr) {
        auto pIt = particles->find(hit.particleId());
        if (pIt != particles->end()) {
          pid = static_cast<std::uint64_t>(
              std::distance(particles->begin(), pIt));
        }
      }
      pidVV->UnsafeAppend(pid);
    }
  }

  auto finish = [](arrow::ListBuilder& b) {
    std::shared_ptr<arrow::Array> out;
    check(b.Finish(&out), "finish list");
    return out;
  };

  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      finish(xList),    finish(yList),   finish(zList),    finish(detList),
      finish(volList),  finish(layList), finish(surfList), finish(pidList),
      finish(txList),   finish(tyList),  finish(tzList),   finish(timeList),
  };

  auto table =
      arrow::Table::Make(ActsPlugins::ArrowUtil::simHitSchema(), arrays);
  m_outputTable(ctx, ActsPlugins::ArrowUtil::ArrowTable{std::move(table)});

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
