// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Arrow/ArrowTrackOutputConverter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
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

ArrowTrackOutputConverter::ArrowTrackOutputConverter(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : ArrowOutputConverter("ArrowTrackOutputConverter", std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing tracks input collection");
  }
  if (m_cfg.outputTable.empty()) {
    throw std::invalid_argument("Missing output table name");
  }
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputTrackParticleMatching.maybeInitialize(
      m_cfg.inputTrackParticleMatching);
  m_inputParticles.maybeInitialize(m_cfg.inputParticles);
  m_inputMeasurementSimHitsMap.maybeInitialize(
      m_cfg.inputMeasurementSimHitsMap);
  m_outputTable.initialize(m_cfg.outputTable);
}

std::vector<std::string> ArrowTrackOutputConverter::collections() const {
  return {m_cfg.outputTable};
}

ProcessCode ArrowTrackOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  const ConstTrackContainer& tracks = m_inputTracks(ctx);
  const TrackParticleMatching* matching =
      m_inputTrackParticleMatching.isInitialized()
          ? &m_inputTrackParticleMatching(ctx)
          : nullptr;
  const SimParticleContainer* particles =
      m_inputParticles.isInitialized() ? &m_inputParticles(ctx) : nullptr;
  const MeasurementSimHitsMap* measToSimHits =
      m_inputMeasurementSimHitsMap.isInitialized()
          ? &m_inputMeasurementSimHitsMap(ctx)
          : nullptr;

  auto* pool = arrow::default_memory_pool();

  arrow::ListBuilder d0List(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder z0List(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder phiList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder thetaList(pool,
                               std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder qopList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder tList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder majIdList(pool,
                               std::make_shared<arrow::UInt64Builder>(pool));
  arrow::ListBuilder hitIdsList(
      pool, std::make_shared<arrow::ListBuilder>(
                pool, std::make_shared<arrow::UInt32Builder>(pool)));
  arrow::ListBuilder trackIdList(pool,
                                 std::make_shared<arrow::UInt16Builder>(pool));

  check(d0List.Append(), "open d0 list");
  check(z0List.Append(), "open z0 list");
  check(phiList.Append(), "open phi list");
  check(thetaList.Append(), "open theta list");
  check(qopList.Append(), "open qop list");
  if (m_cfg.writeTime) {
    check(tList.Append(), "open t list");
  } else {
    // Whole-list null: encodes "no time data" with one definition-level
    // entry instead of N per-track inner nulls.
    check(tList.AppendNull(), "append null t list");
  }
  check(majIdList.Append(), "open majority_particle_id list");
  check(hitIdsList.Append(), "open hit_ids outer list");
  check(trackIdList.Append(), "open track_id list");

  auto* d0V = static_cast<arrow::FloatBuilder*>(d0List.value_builder());
  auto* z0V = static_cast<arrow::FloatBuilder*>(z0List.value_builder());
  auto* phiV = static_cast<arrow::FloatBuilder*>(phiList.value_builder());
  auto* thetaV = static_cast<arrow::FloatBuilder*>(thetaList.value_builder());
  auto* qopV = static_cast<arrow::FloatBuilder*>(qopList.value_builder());
  auto* tV = static_cast<arrow::FloatBuilder*>(tList.value_builder());
  auto* majIdV = static_cast<arrow::UInt64Builder*>(majIdList.value_builder());
  auto* hitIdsInner =
      static_cast<arrow::ListBuilder*>(hitIdsList.value_builder());
  auto* hitIdsV =
      static_cast<arrow::UInt32Builder*>(hitIdsInner->value_builder());
  auto* trackIdV =
      static_cast<arrow::UInt16Builder*>(trackIdList.value_builder());

  const auto n = tracks.size();
  check(d0V->Reserve(n), "reserve d0");
  check(z0V->Reserve(n), "reserve z0");
  check(phiV->Reserve(n), "reserve phi");
  check(thetaV->Reserve(n), "reserve theta");
  check(qopV->Reserve(n), "reserve qop");
  if (m_cfg.writeTime) {
    check(tV->Reserve(n), "reserve t");
  }
  check(majIdV->Reserve(n), "reserve majority_particle_id");
  check(trackIdV->Reserve(n), "reserve track_id");

  // Sentinel for "no matched particle": an out-of-range row index. Real
  // indices are < particles->size(), so this can never collide.
  // @TODO: Turn into explicit optionals
  constexpr std::uint64_t kUnmatched =
      std::numeric_limits<std::uint64_t>::max();

  for (const auto& track : tracks) {
    if (track.hasReferenceSurface()) {
      const auto& p = track.parameters();
      d0V->UnsafeAppend(static_cast<float>(p[Acts::eBoundLoc0]));
      z0V->UnsafeAppend(static_cast<float>(p[Acts::eBoundLoc1]));
      phiV->UnsafeAppend(static_cast<float>(p[Acts::eBoundPhi]));
      thetaV->UnsafeAppend(static_cast<float>(p[Acts::eBoundTheta]));
      qopV->UnsafeAppend(static_cast<float>(p[Acts::eBoundQOverP]));
      if (m_cfg.writeTime) {
        tV->UnsafeAppend(static_cast<float>(p[Acts::eBoundTime]));
      }
    } else {
      d0V->UnsafeAppendNull();
      z0V->UnsafeAppendNull();
      phiV->UnsafeAppendNull();
      thetaV->UnsafeAppendNull();
      qopV->UnsafeAppendNull();
      if (m_cfg.writeTime) {
        tV->UnsafeAppendNull();
      }
    }

    // Resolve the matched truth barcode to its row index in the configured
    // particle container so it joins against the particle table's
    // `particle_id` column.
    std::uint64_t majId = kUnmatched;
    if (matching != nullptr && particles != nullptr) {
      auto it = matching->find(track.index());
      if (it != matching->end() && it->second.particle.has_value()) {
        const auto& bc = it->second.particle.value();
        auto pIt = particles->find(bc);
        if (pIt != particles->end()) {
          majId = static_cast<std::uint64_t>(
              std::distance(particles->begin(), pIt));
        }
      }
    }
    majIdV->UnsafeAppend(majId);

    check(hitIdsInner->Append(), "open hit_ids inner list");
    // Without a measurement→sim-hit map we leave the inner list empty rather
    // than silently emit measurement indices, which would alias mismatched
    // ids into a downstream sim-hit table.
    if (measToSimHits != nullptr) {
      // `trackStatesReversed()` walks outermost→innermost (the only direct
      // iteration the MultiTrajectory proxy offers); buffer and reverse so
      // the per-track list comes out inner→outer along the trajectory, which
      // is what downstream consumers expect.
      std::vector<std::uint32_t> hitIds;
      for (const auto& state : track.trackStatesReversed()) {
        if (!state.hasUncalibratedSourceLink()) {
          continue;
        }
        const auto sl =
            state.getUncalibratedSourceLink().template get<IndexSourceLink>();
        const auto measIdx = static_cast<Index>(sl.index());
        // One measurement may map to multiple sim hits (clustering merged
        // them); flatten them into the per-track list.
        auto range = measToSimHits->equal_range(measIdx);
        for (auto it = range.first; it != range.second; ++it) {
          hitIds.push_back(static_cast<std::uint32_t>(it->second));
        }
      }
      for (auto rit = hitIds.rbegin(); rit != hitIds.rend(); ++rit) {
        check(hitIdsV->Append(*rit), "append hit_id");
      }
    }

    trackIdV->UnsafeAppend(static_cast<std::uint16_t>(track.index()));
  }

  auto finish = [](arrow::ListBuilder& b) {
    std::shared_ptr<arrow::Array> out;
    check(b.Finish(&out), "finish list");
    return out;
  };

  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      finish(d0List),     finish(z0List),      finish(phiList),
      finish(thetaList),  finish(qopList),     finish(majIdList),
      finish(hitIdsList), finish(trackIdList), finish(tList),
  };

  auto table =
      arrow::Table::Make(ActsPlugins::ArrowUtil::trackSchema(), arrays);
  m_outputTable(ctx, ActsPlugins::ArrowUtil::ArrowTable{std::move(table)});

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
