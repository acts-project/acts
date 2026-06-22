// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Arrow/ArrowMeasurementOutputConverter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

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

ArrowMeasurementOutputConverter::ArrowMeasurementOutputConverter(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : ArrowOutputConverter("ArrowMeasurementOutputConverter",
                           std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
  }
  // Truth is optional (e.g. data, or reco-only conversions), but the sim-hit
  // container and the sim-hit -> measurement map only make sense as a pair.
  if (m_cfg.inputSimHits.empty() != m_cfg.inputSimHitMeasurementsMap.empty()) {
    throw std::invalid_argument(
        "inputSimHits and inputSimHitMeasurementsMap must be set together");
  }
  if (m_cfg.outputTable.empty()) {
    throw std::invalid_argument("Missing output table name");
  }
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("trackingGeometry must be set");
  }
  if (!m_cfg.detectorResolver) {
    throw std::invalid_argument("detectorResolver must be set");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
  m_inputSimHits.maybeInitialize(m_cfg.inputSimHits);
  m_inputParticles.maybeInitialize(m_cfg.inputParticles);
  m_inputSimHitMeasurementsMap.maybeInitialize(m_cfg.inputSimHitMeasurementsMap);
  m_outputTable.initialize(m_cfg.outputTable);
}

std::vector<std::string> ArrowMeasurementOutputConverter::collections() const {
  return {m_cfg.outputTable};
}

ProcessCode ArrowMeasurementOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  const MeasurementContainer& measurements = m_inputMeasurements(ctx);
  const ClusterContainer* clusters =
      m_inputClusters.isInitialized() ? &m_inputClusters(ctx) : nullptr;
  const SimHitContainer* simHits =
      m_inputSimHits.isInitialized() ? &m_inputSimHits(ctx) : nullptr;
  const SimParticleContainer* particles =
      m_inputParticles.isInitialized() ? &m_inputParticles(ctx) : nullptr;
  const SimHitMeasurementsMap* simHitMeasMap =
      m_inputSimHitMeasurementsMap.isInitialized()
          ? &m_inputSimHitMeasurementsMap(ctx)
          : nullptr;

  if (clusters != nullptr && clusters->size() != measurements.size()) {
    throw std::runtime_error(
        "ArrowMeasurementOutputConverter: clusters (" +
        std::to_string(clusters->size()) + ") and measurements (" +
        std::to_string(measurements.size()) +
        ") are not parallel containers");
  }

  // Invert the sim-hit -> measurement map into measurement -> sim-hits.
  // Sim-hits referencing no measurement stay in the truth table; here they
  // only feed a diagnostic.
  std::unordered_multimap<Index, SimHitIndex> measToSimHits;
  if (simHitMeasMap != nullptr) {
    measToSimHits.reserve(simHitMeasMap->size());
    std::vector<bool> referenced(simHits->size(), false);
    for (const auto& [simHitIdx, measIdx] : *simHitMeasMap) {
      measToSimHits.emplace(measIdx, simHitIdx);
      if (simHitIdx < referenced.size()) {
        referenced[simHitIdx] = true;
      }
    }
    std::size_t droppedSimHits = 0;
    for (bool seen : referenced) {
      if (!seen) {
        ++droppedSimHits;
      }
    }
    const double droppedFraction =
        simHits->empty() ? 0.0
                         : static_cast<double>(droppedSimHits) /
                               static_cast<double>(simHits->size());
    if (droppedFraction > m_cfg.maxUnmatchedSimHitFraction) {
      ACTS_WARNING("ArrowMeasurementOutputConverter: "
                   << droppedSimHits << " / " << simHits->size()
                   << " sim-hits (" << droppedFraction * 100.0
                   << "%) contribute to no measurement (links absent from this "
                      "table; the hits remain in the truth table)");
    }
  }

  auto* pool = arrow::default_memory_pool();

  // Flat per-measurement columns. (Arrow builders are neither copyable nor
  // movable, so each is constructed in place.)
  arrow::ListBuilder loc0List(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder loc1List(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vl0List(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vl1List(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder timeList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vtList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder subList(pool, std::make_shared<arrow::UInt8Builder>(pool));
  arrow::ListBuilder xList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder yList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder zList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder detList(pool, std::make_shared<arrow::UInt8Builder>(pool));
  arrow::ListBuilder volList(pool, std::make_shared<arrow::UInt8Builder>(pool));
  arrow::ListBuilder layList(pool,
                             std::make_shared<arrow::UInt16Builder>(pool));
  arrow::ListBuilder surfList(pool,
                              std::make_shared<arrow::UInt32Builder>(pool));
  arrow::ListBuilder sz0List(pool,
                             std::make_shared<arrow::UInt16Builder>(pool));
  arrow::ListBuilder sz1List(pool,
                             std::make_shared<arrow::UInt16Builder>(pool));
  arrow::ListBuilder nchList(pool,
                             std::make_shared<arrow::UInt16Builder>(pool));
  arrow::ListBuilder sumList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder letaList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder lphiList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder getaList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder gphiList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder eangList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pangList(pool, std::make_shared<arrow::FloatBuilder>(pool));

  // Nested truth-link columns: outer = per measurement, inner = per
  // contributing sim-hit.
  auto pidInner = std::make_shared<arrow::ListBuilder>(
      pool, std::make_shared<arrow::UInt64Builder>(pool));
  auto shidInner = std::make_shared<arrow::ListBuilder>(
      pool, std::make_shared<arrow::UInt32Builder>(pool));
  arrow::ListBuilder pidList(pool, pidInner);
  arrow::ListBuilder shidList(pool, shidInner);

  std::vector<arrow::ListBuilder*> all = {
      &loc0List, &loc1List, &vl0List, &vl1List, &timeList, &vtList,
      &subList,  &xList,    &yList,   &zList,   &detList,  &volList,
      &layList,  &surfList, &sz0List, &sz1List, &nchList,  &sumList,
      &letaList, &lphiList, &getaList, &gphiList, &eangList, &pangList,
      &pidList,  &shidList};
  for (auto* b : all) {
    check(b->Append(), "open per-event list");
  }

  auto fv = [](arrow::ListBuilder& b) {
    return static_cast<arrow::FloatBuilder*>(b.value_builder());
  };
  auto* loc0V = fv(loc0List);
  auto* loc1V = fv(loc1List);
  auto* vl0V = fv(vl0List);
  auto* vl1V = fv(vl1List);
  auto* timeV = fv(timeList);
  auto* vtV = fv(vtList);
  auto* subV = static_cast<arrow::UInt8Builder*>(subList.value_builder());
  auto* xV = fv(xList);
  auto* yV = fv(yList);
  auto* zV = fv(zList);
  auto* detV = static_cast<arrow::UInt8Builder*>(detList.value_builder());
  auto* volV = static_cast<arrow::UInt8Builder*>(volList.value_builder());
  auto* layV = static_cast<arrow::UInt16Builder*>(layList.value_builder());
  auto* surfV = static_cast<arrow::UInt32Builder*>(surfList.value_builder());
  auto* sz0V = static_cast<arrow::UInt16Builder*>(sz0List.value_builder());
  auto* sz1V = static_cast<arrow::UInt16Builder*>(sz1List.value_builder());
  auto* nchV = static_cast<arrow::UInt16Builder*>(nchList.value_builder());
  auto* sumV = fv(sumList);
  auto* letaV = fv(letaList);
  auto* lphiV = fv(lphiList);
  auto* getaV = fv(getaList);
  auto* gphiV = fv(gphiList);
  auto* eangV = fv(eangList);
  auto* pangV = fv(pangList);

  auto* pidVL = static_cast<arrow::ListBuilder*>(pidList.value_builder());
  auto* shidVL = static_cast<arrow::ListBuilder*>(shidList.value_builder());
  auto* pidVV = static_cast<arrow::UInt64Builder*>(pidVL->value_builder());
  auto* shidVV = static_cast<arrow::UInt32Builder*>(shidVL->value_builder());

  constexpr std::uint64_t kUnmatched =
      std::numeric_limits<std::uint64_t>::max();
  constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();

  auto clampU16 = [](std::size_t v) {
    return static_cast<std::uint16_t>(
        std::min<std::size_t>(v, std::numeric_limits<std::uint16_t>::max()));
  };

  const std::size_t nMeas = measurements.size();
  // Reserve every flat value builder before the UnsafeAppend loop — without
  // this, UnsafeAppend writes past the buffer (heap corruption, not a clean
  // failure).
  for (auto* fb : {loc0V, loc1V, vl0V, vl1V, timeV, vtV, xV, yV, zV, sumV,
                   letaV, lphiV, getaV, gphiV, eangV, pangV}) {
    check(fb->Reserve(nMeas), "reserve float column");
  }
  check(subV->Reserve(nMeas), "reserve subspace");
  check(detV->Reserve(nMeas), "reserve detector");
  check(volV->Reserve(nMeas), "reserve volume_id");
  for (auto* ub : {layV, sz0V, sz1V, nchV}) {
    check(ub->Reserve(nMeas), "reserve uint16 column");
  }
  check(surfV->Reserve(nMeas), "reserve surface_id");

  for (Index measIdx = 0; measIdx < nMeas; ++measIdx) {
    const auto meas = measurements.getMeasurement(measIdx);
    const auto gid = meas.geometryId();

    // Local parameters/variances, expanded to the stable full bound layout.
    // The subspace bitmask records which components were actually measured.
    const auto full = meas.fullParameters();
    const auto cov = meas.fullCovariance();
    loc0V->UnsafeAppend(static_cast<float>(full[Acts::eBoundLoc0] /
                                           Acts::UnitConstants::mm));
    loc1V->UnsafeAppend(static_cast<float>(full[Acts::eBoundLoc1] /
                                           Acts::UnitConstants::mm));
    vl0V->UnsafeAppend(static_cast<float>(
        cov(Acts::eBoundLoc0, Acts::eBoundLoc0) /
        (Acts::UnitConstants::mm * Acts::UnitConstants::mm)));
    vl1V->UnsafeAppend(static_cast<float>(
        cov(Acts::eBoundLoc1, Acts::eBoundLoc1) /
        (Acts::UnitConstants::mm * Acts::UnitConstants::mm)));
    timeV->UnsafeAppend(static_cast<float>(full[Acts::eBoundTime] /
                                           Acts::UnitConstants::mm));
    vtV->UnsafeAppend(static_cast<float>(
        cov(Acts::eBoundTime, Acts::eBoundTime) /
        (Acts::UnitConstants::mm * Acts::UnitConstants::mm)));

    std::uint8_t subspace = 0;
    if (meas.contains(Acts::eBoundLoc0)) {
      subspace |= 1u;
    }
    if (meas.contains(Acts::eBoundLoc1)) {
      subspace |= 2u;
    }
    if (meas.contains(Acts::eBoundTime)) {
      subspace |= 4u;
    }
    subV->UnsafeAppend(subspace);

    // Reco (digitized) global position: project the measurement's bound
    // parameters through its surface. NaN for non-regular surfaces.
    float gx = kNaN;
    float gy = kNaN;
    float gz = kNaN;
    const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(gid);
    if (surface == nullptr) {
      throw std::runtime_error(
          "ArrowMeasurementOutputConverter: surface not found for geometry "
          "id " +
          std::to_string(gid.value()));
    }
    if (const auto* regular =
            dynamic_cast<const Acts::RegularSurface*>(surface)) {
      Acts::Vector2 loc{full[Acts::eBoundLoc0], full[Acts::eBoundLoc1]};
      const Acts::Vector3 global = regular->localToGlobal(ctx.geoContext, loc);
      gx = static_cast<float>(global.x() / Acts::UnitConstants::mm);
      gy = static_cast<float>(global.y() / Acts::UnitConstants::mm);
      gz = static_cast<float>(global.z() / Acts::UnitConstants::mm);
    }
    xV->UnsafeAppend(gx);
    yV->UnsafeAppend(gy);
    zV->UnsafeAppend(gz);

    detV->UnsafeAppend(m_cfg.detectorResolver(gid));
    volV->UnsafeAppend(static_cast<std::uint8_t>(gid.volume()));
    layV->UnsafeAppend(static_cast<std::uint16_t>(gid.layer()));
    surfV->UnsafeAppend(static_cast<std::uint32_t>(gid.sensitive()));

    // Cluster-shape features (parallel container, checked above). Without a
    // clusters input the shape columns are zeros - the schema stays stable.
    static const Cluster kEmptyCluster{};
    const Cluster& clu =
        clusters != nullptr ? (*clusters)[measIdx] : kEmptyCluster;
    sz0V->UnsafeAppend(clampU16(clu.sizeLoc0));
    sz1V->UnsafeAppend(clampU16(clu.sizeLoc1));
    nchV->UnsafeAppend(clampU16(clu.channels.size()));
    sumV->UnsafeAppend(static_cast<float>(clu.sumActivations()));
    letaV->UnsafeAppend(clu.localEta);
    lphiV->UnsafeAppend(clu.localPhi);
    getaV->UnsafeAppend(clu.globalEta);
    gphiV->UnsafeAppend(clu.globalPhi);
    eangV->UnsafeAppend(clu.etaAngle);
    pangV->UnsafeAppend(clu.phiAngle);

    // Truth links: contributing sim-hits (row indices into the truth table)
    // and their particles (row indices into the particle table).
    check(pidVL->Append(), "open per-measurement particle_ids list");
    check(shidVL->Append(), "open per-measurement simhit_ids list");
    auto range = measToSimHits.equal_range(measIdx);
    const auto nContribs =
        static_cast<std::int64_t>(std::distance(range.first, range.second));
    check(pidVV->Reserve(nContribs), "reserve particle_ids contribs");
    check(shidVV->Reserve(nContribs), "reserve simhit_ids contribs");
    for (auto it = range.first; it != range.second; ++it) {
      const SimHitIndex simHitIdx = it->second;
      shidVV->UnsafeAppend(static_cast<std::uint32_t>(simHitIdx));

      std::uint64_t pid = kUnmatched;
      if (particles != nullptr && simHits != nullptr) {
        const auto& hit = *simHits->nth(simHitIdx);
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
      finish(loc0List), finish(loc1List), finish(vl0List),  finish(vl1List),
      finish(timeList), finish(vtList),   finish(subList),  finish(xList),
      finish(yList),    finish(zList),    finish(detList),  finish(volList),
      finish(layList),  finish(surfList), finish(sz0List),  finish(sz1List),
      finish(nchList),  finish(sumList),  finish(letaList), finish(lphiList),
      finish(getaList), finish(gphiList), finish(eangList), finish(pangList),
      finish(pidList),  finish(shidList),
  };

  auto table =
      arrow::Table::Make(ActsPlugins::ArrowUtil::measurementSchema(), arrays);
  m_outputTable(ctx, ActsPlugins::ArrowUtil::ArrowTable{std::move(table)});

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
