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

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>

#include <arrow/api.h>

namespace ActsExamples {

namespace {

std::shared_ptr<arrow::Schema> simHitSchema() {
  return arrow::schema({
      arrow::field("x", arrow::list(arrow::float32()), false),
      arrow::field("y", arrow::list(arrow::float32()), false),
      arrow::field("z", arrow::list(arrow::float32()), false),
      arrow::field("true_x", arrow::list(arrow::float32()), false),
      arrow::field("true_y", arrow::list(arrow::float32()), false),
      arrow::field("true_z", arrow::list(arrow::float32()), false),
      arrow::field("time", arrow::list(arrow::float32()), false),
      arrow::field("particle_id", arrow::list(arrow::uint64()), false),
      arrow::field("detector", arrow::list(arrow::uint8()), false),
      arrow::field("volume_id", arrow::list(arrow::uint8()), false),
      arrow::field("layer_id", arrow::list(arrow::uint16()), false),
      arrow::field("surface_id", arrow::list(arrow::uint32()), false),
  });
}

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

  // The digitized x,y,z columns require the trio (measurements, map, geometry)
  // — partial wiring would silently produce stale or NaN positions, so reject
  // it up front.
  const bool hasMeas = !m_cfg.inputMeasurements.empty();
  const bool hasMap = !m_cfg.inputSimHitMeasurementsMap.empty();
  const bool hasGeo = m_cfg.trackingGeometry != nullptr;
  if ((hasMeas || hasMap || hasGeo) && !(hasMeas && hasMap && hasGeo)) {
    throw std::invalid_argument(
        "ArrowSimHitOutputConverter: inputMeasurements, "
        "inputSimHitMeasurementsMap, and trackingGeometry must all be set or "
        "all be unset");
  }

  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputParticles.maybeInitialize(m_cfg.inputParticles);
  m_inputMeasurements.maybeInitialize(m_cfg.inputMeasurements);
  m_inputSimHitMeasurementsMap.maybeInitialize(
      m_cfg.inputSimHitMeasurementsMap);
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
  const MeasurementContainer* measurements =
      m_inputMeasurements.isInitialized() ? &m_inputMeasurements(ctx) : nullptr;
  const SimHitMeasurementsMap* simHitMeasMap =
      m_inputSimHitMeasurementsMap.isInitialized()
          ? &m_inputSimHitMeasurementsMap(ctx)
          : nullptr;

  auto* pool = arrow::default_memory_pool();

  arrow::ListBuilder xList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder yList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder zList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder txList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder tyList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder tzList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder timeList(pool,
                              std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pidList(pool,
                             std::make_shared<arrow::UInt64Builder>(pool));
  arrow::ListBuilder detList(pool, std::make_shared<arrow::UInt8Builder>(pool));
  arrow::ListBuilder volList(pool, std::make_shared<arrow::UInt8Builder>(pool));
  arrow::ListBuilder layList(pool,
                             std::make_shared<arrow::UInt16Builder>(pool));
  arrow::ListBuilder surfList(pool,
                              std::make_shared<arrow::UInt32Builder>(pool));

  check(xList.Append(), "open x list");
  check(yList.Append(), "open y list");
  check(zList.Append(), "open z list");
  check(txList.Append(), "open true_x list");
  check(tyList.Append(), "open true_y list");
  check(tzList.Append(), "open true_z list");
  check(timeList.Append(), "open time list");
  check(pidList.Append(), "open particle_id list");
  check(detList.Append(), "open detector list");
  check(volList.Append(), "open volume_id list");
  check(layList.Append(), "open layer_id list");
  check(surfList.Append(), "open surface_id list");

  // TODO: Keep typed child builder handles when constructing the list builders
  // instead of recovering them through value_builder() and static_cast.
  auto* xV = static_cast<arrow::FloatBuilder*>(xList.value_builder());
  auto* yV = static_cast<arrow::FloatBuilder*>(yList.value_builder());
  auto* zV = static_cast<arrow::FloatBuilder*>(zList.value_builder());
  auto* txV = static_cast<arrow::FloatBuilder*>(txList.value_builder());
  auto* tyV = static_cast<arrow::FloatBuilder*>(tyList.value_builder());
  auto* tzV = static_cast<arrow::FloatBuilder*>(tzList.value_builder());
  auto* timeV = static_cast<arrow::FloatBuilder*>(timeList.value_builder());
  auto* pidV = static_cast<arrow::UInt64Builder*>(pidList.value_builder());
  auto* detV = static_cast<arrow::UInt8Builder*>(detList.value_builder());
  auto* volV = static_cast<arrow::UInt8Builder*>(volList.value_builder());
  auto* layV = static_cast<arrow::UInt16Builder*>(layList.value_builder());
  auto* surfV = static_cast<arrow::UInt32Builder*>(surfList.value_builder());

  const auto n = simHits.size();
  check(xV->Reserve(n), "reserve x");
  check(yV->Reserve(n), "reserve y");
  check(zV->Reserve(n), "reserve z");
  check(txV->Reserve(n), "reserve true_x");
  check(tyV->Reserve(n), "reserve true_y");
  check(tzV->Reserve(n), "reserve true_z");
  check(timeV->Reserve(n), "reserve time");
  check(pidV->Reserve(n), "reserve particle_id");
  check(detV->Reserve(n), "reserve detector");
  check(volV->Reserve(n), "reserve volume_id");
  check(layV->Reserve(n), "reserve layer_id");
  check(surfV->Reserve(n), "reserve surface_id");

  // Sentinel matches the convention used by ArrowTrackOutputConverter for
  // unmatched rows.
  // @TODO: Turn into explicit optionals?
  constexpr std::uint64_t kUnmatched =
      std::numeric_limits<std::uint64_t>::max();
  constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();

  SimHitIndex hitIdx = 0;
  for (const auto& hit : simHits) {
    const auto& pos = hit.fourPosition();
    const float tx = static_cast<float>(pos.x() / Acts::UnitConstants::mm);
    const float ty = static_cast<float>(pos.y() / Acts::UnitConstants::mm);
    const float tz = static_cast<float>(pos.z() / Acts::UnitConstants::mm);
    const float t = static_cast<float>(pos.w() / Acts::UnitConstants::mm);

    txV->UnsafeAppend(tx);
    tyV->UnsafeAppend(ty);
    tzV->UnsafeAppend(tz);
    timeV->UnsafeAppend(t);

    // Digitized global position: project the first matched measurement's
    // bound parameters through the surface. Multiple measurements per hit
    // would only happen if a hit migrated across modules during clustering;
    // we take the first deterministically and leave the rest for a future
    // "merged hits" extension.
    float gx = kNaN;
    float gy = kNaN;
    float gz = kNaN;
    if (simHitMeasMap != nullptr && measurements != nullptr) {
      auto range = simHitMeasMap->equal_range(hitIdx);
      if (range.first != range.second) {
        const Index measIdx = range.first->second;
        const auto meas = measurements->getMeasurement(measIdx);
        const Acts::Surface* surface =
            m_cfg.trackingGeometry->findSurface(meas.geometryId());
        if (surface == nullptr) {
          throw std::runtime_error(
              "ArrowSimHitOutputConverter: surface not found for geometry id " +
              std::to_string(meas.geometryId().value()));
        }
        const auto* regular =
            dynamic_cast<const Acts::RegularSurface*>(surface);
        if (regular != nullptr) {
          // Expand the (possibly 1D) measurement to a full bound vector so
          // local positions are read from a stable layout regardless of which
          // subspace coordinates were measured.
          const auto full = meas.fullParameters();
          Acts::Vector2 loc{full[Acts::eBoundLoc0], full[Acts::eBoundLoc1]};
          const Acts::Vector3 global =
              regular->localToGlobal(ctx.geoContext, loc);
          gx = static_cast<float>(global.x() / Acts::UnitConstants::mm);
          gy = static_cast<float>(global.y() / Acts::UnitConstants::mm);
          gz = static_cast<float>(global.z() / Acts::UnitConstants::mm);
        }
      }
    }
    xV->UnsafeAppend(gx);
    yV->UnsafeAppend(gy);
    zV->UnsafeAppend(gz);

    std::uint64_t pid = kUnmatched;
    if (particles != nullptr) {
      auto pIt = particles->find(hit.particleId());
      if (pIt != particles->end()) {
        pid =
            static_cast<std::uint64_t>(std::distance(particles->begin(), pIt));
      }
    }
    pidV->UnsafeAppend(pid);

    const auto gid = hit.geometryId();
    volV->UnsafeAppend(static_cast<std::uint8_t>(gid.volume()));
    layV->UnsafeAppend(static_cast<std::uint16_t>(gid.layer()));
    surfV->UnsafeAppend(static_cast<std::uint32_t>(gid.sensitive()));
    // The default `detectorResolver` reads the geometry id's
    // `extra` byte, which we rely on geometry construction to stamp with a
    // per-surface subsystem id. By default, every hit gets `extra() == 0`
    // unless the user supplies a custom resolver.
    detV->UnsafeAppend(m_cfg.detectorResolver(gid));

    ++hitIdx;
  }

  auto finish = [](arrow::ListBuilder& b) {
    std::shared_ptr<arrow::Array> out;
    check(b.Finish(&out), "finish list");
    return out;
  };

  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      finish(xList),   finish(yList),   finish(zList),    finish(txList),
      finish(tyList),  finish(tzList),  finish(timeList), finish(pidList),
      finish(detList), finish(volList), finish(layList),  finish(surfList),
  };

  auto table = arrow::Table::Make(simHitSchema(), arrays);
  m_outputTable(ctx, std::move(table));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
