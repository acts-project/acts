// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/TwoTrackVertexFinder.hpp"

#include "Acts/MagneticField/MagneticFieldQuery.hpp"
#include "Acts/Vertexing/IVertexFitter.hpp"
#include "Acts/Vertexing/VertexingError.hpp"
#include "Acts/Vertexing/detail/TwoTrackVertexSeed.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>

#include <Eigen/Cholesky>

namespace Acts {
namespace {
struct PairFinderState {
  MagneticFieldQuery field;
  std::vector<InputTrack> removed;
};
}  // namespace

TwoTrackVertexFinder::TwoTrackVertexFinder(Config cfg) : m_cfg(std::move(cfg)) {
  if (!m_cfg.vertexFitter || !m_cfg.field ||
      !m_cfg.extractParameters.connected() ||
      !std::isfinite(m_cfg.maxSeedSeparation) || m_cfg.maxSeedSeparation < 0. ||
      !std::isfinite(m_cfg.maxChi2) || m_cfg.maxChi2 < 0.) {
    throw std::invalid_argument("TwoTrackVertexFinder: invalid configuration");
  }
}

IVertexFinder::State TwoTrackVertexFinder::makeState(
    const MagneticFieldContext& mctx) const {
  return IVertexFinder::State{
      PairFinderState{MagneticFieldQuery(*m_cfg.field, mctx), {}}};
}

void TwoTrackVertexFinder::setTracksToRemove(
    IVertexFinder::State& state, const std::vector<InputTrack>& tracks) const {
  auto& removed = state.as<PairFinderState>().removed;
  removed.insert(removed.end(), tracks.begin(), tracks.end());
}

Result<std::vector<Vertex>> TwoTrackVertexFinder::find(
    const std::vector<InputTrack>& tracks, const VertexingOptions& options,
    IVertexFinder::State& state) const {
  auto& data = state.as<PairFinderState>();
  if (&options.magFieldContext.get() != &data.field.context()) {
    return VertexingError::InvalidInput;
  }
  if (options.useConstraintInFit) {
    // Billoir consumes the full spacetime covariance, including time.
    const auto& covariance = options.constraint.fullCovariance();
    if (!options.constraint.fullPosition().allFinite() ||
        !covariance.allFinite() ||
        !covariance.isApprox(covariance.transpose()) ||
        covariance.llt().info() != Eigen::Success) {
      return VertexingError::InvalidInput;
    }
  }
  const std::optional<Vertex> constraint =
      options.useConstraintInFit ? std::optional<Vertex>(options.constraint)
                                 : std::nullopt;
  std::vector<InputTrack> active;
  std::vector<BoundTrackParameters> parameters;
  for (const auto& track : tracks) {
    if (std::ranges::find(data.removed, track) == data.removed.end() &&
        std::ranges::find(active, track) == active.end()) {
      active.push_back(track);
      parameters.push_back(m_cfg.extractParameters(track));
    }
  }
  std::vector<Vertex> vertices;
  for (std::size_t i = 0; i < active.size(); ++i) {
    for (std::size_t j = i + 1; j < active.size(); ++j) {
      const auto& first = parameters[i];
      const auto& second = parameters[j];
      if (m_cfg.oppositeChargeOnly && first.charge() * second.charge() >= 0.) {
        continue;
      }
      const Vector3 midpoint = (first.position(options.geoContext) +
                                second.position(options.geoContext)) /
                               2.;
      auto field = data.field.getField(midpoint);
      if (!field.ok()) {
        return field.error();
      }
      if (!field->allFinite() ||
          field->head<2>().norm() > 1e-6 * std::abs(field->z())) {
        return VertexingError::InvalidInput;
      }
      auto seed = detail::estimateTwoTrackVertexSeed(
          first, second, options.geoContext, field->z());
      if (!seed.ok()) {
        if (seed.error() == VertexingError::SeedingError) {
          continue;
        }
        return seed.error();
      }
      if (seed->separation > m_cfg.maxSeedSeparation) {
        continue;
      }
      Vector4 position;
      position.head<3>() = seed->position;
      // Spatial prototype: the fitter determines time from the track
      // parameters.
      position[eTime] = (first.time() + second.time()) / 2.;
      const std::array pair{active[i], active[j]};
      auto vertex =
          m_cfg.vertexFitter->fit(VertexFitInput{pair, position, constraint},
                                  options.geoContext, data.field);
      if (!vertex.ok()) {
        return vertex.error();
      }
      if (vertex->fullPosition().allFinite() &&
          vertex->fullCovariance().allFinite() &&
          std::isfinite(vertex->fitQuality().first) &&
          vertex->fitQuality().first >= 0. &&
          vertex->fitQuality().first <= m_cfg.maxChi2) {
        vertices.push_back(std::move(*vertex));
      }
    }
  }
  return vertices;
}
}  // namespace Acts
