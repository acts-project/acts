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

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>

#include <arrow/api.h>

namespace ActsExamples {

namespace {

/// Nested layout: one row per event, each field a @c list<T> whose single
/// list element at row @c N holds all tracks of event @c N. The @c hit_ids
/// column is a list-of-list so each track keeps its own vector of
/// measurement indices.
std::shared_ptr<arrow::Schema> trackSchema() {
  return arrow::schema({
      arrow::field("d0", arrow::list(arrow::float32()), false),
      arrow::field("z0", arrow::list(arrow::float32()), false),
      arrow::field("phi", arrow::list(arrow::float32()), false),
      arrow::field("theta", arrow::list(arrow::float32()), false),
      arrow::field("qop", arrow::list(arrow::float32()), false),
      arrow::field("t", arrow::list(arrow::float32()), false),
      arrow::field("majority_particle_id", arrow::list(arrow::uint64()), false),
      arrow::field("hit_ids", arrow::list(arrow::list(arrow::uint32())), false),
      arrow::field("track_id", arrow::list(arrow::uint16()), false),
  });
}

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
  m_inputTrackParticleMatching.maybeInitialize(m_cfg.inputTrackParticleMatching);
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
  arrow::ListBuilder trackIdList(
      pool, std::make_shared<arrow::UInt16Builder>(pool));

  check(d0List.Append(), "open d0 list");
  check(z0List.Append(), "open z0 list");
  check(phiList.Append(), "open phi list");
  check(thetaList.Append(), "open theta list");
  check(qopList.Append(), "open qop list");
  check(tList.Append(), "open t list");
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
  auto* hitIdsV = static_cast<arrow::UInt32Builder*>(hitIdsInner->value_builder());
  auto* trackIdV =
      static_cast<arrow::UInt16Builder*>(trackIdList.value_builder());

  const auto n = tracks.size();
  check(d0V->Reserve(n), "reserve d0");
  check(z0V->Reserve(n), "reserve z0");
  check(phiV->Reserve(n), "reserve phi");
  check(thetaV->Reserve(n), "reserve theta");
  check(qopV->Reserve(n), "reserve qop");
  check(tV->Reserve(n), "reserve t");
  check(majIdV->Reserve(n), "reserve majority_particle_id");
  check(trackIdV->Reserve(n), "reserve track_id");

  constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();

  for (const auto& track : tracks) {
    if (track.hasReferenceSurface()) {
      const auto& p = track.parameters();
      d0V->UnsafeAppend(static_cast<float>(p[Acts::eBoundLoc0]));
      z0V->UnsafeAppend(static_cast<float>(p[Acts::eBoundLoc1]));
      phiV->UnsafeAppend(static_cast<float>(p[Acts::eBoundPhi]));
      thetaV->UnsafeAppend(static_cast<float>(p[Acts::eBoundTheta]));
      qopV->UnsafeAppend(static_cast<float>(p[Acts::eBoundQOverP]));
      tV->UnsafeAppend(static_cast<float>(p[Acts::eBoundTime]));
    } else {
      d0V->UnsafeAppend(kNaN);
      z0V->UnsafeAppend(kNaN);
      phiV->UnsafeAppend(kNaN);
      thetaV->UnsafeAppend(kNaN);
      qopV->UnsafeAppend(kNaN);
      tV->UnsafeAppend(kNaN);
    }

    std::uint64_t majId = 0;
    if (matching != nullptr) {
      auto it = matching->find(track.index());
      if (it != matching->end() && it->second.particle.has_value()) {
        majId = static_cast<std::uint64_t>(
            it->second.particle.value().particle());
      }
    }
    majIdV->UnsafeAppend(majId);

    check(hitIdsInner->Append(), "open hit_ids inner list");
    for (const auto& state : track.trackStatesReversed()) {
      if (!state.hasUncalibratedSourceLink()) {
        continue;
      }
      const auto sl =
          state.getUncalibratedSourceLink().template get<IndexSourceLink>();
      check(hitIdsV->Append(static_cast<std::uint32_t>(sl.index())),
            "append hit_id");
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
      finish(thetaList),  finish(qopList),     finish(tList),
      finish(majIdList),  finish(hitIdsList),  finish(trackIdList),
  };

  auto table = arrow::Table::Make(trackSchema(), arrays);
  m_outputTable(ctx, std::move(table));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
