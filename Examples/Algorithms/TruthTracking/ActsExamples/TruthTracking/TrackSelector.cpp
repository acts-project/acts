// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TrackSelector.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/TruthTracking/TruthVerticesToTracks.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

ActsExamples::TrackSelector::TrackSelector(const Config& cfg,
                                           Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("Selector", level), m_cfg(cfg) {
  if (m_cfg.input.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
  if (m_cfg.output.empty()) {
    throw std::invalid_argument("Missing output collection");
  }
}

ActsExamples::ProcessCode ActsExamples::TrackSelector::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  std::vector<VertexAndTracks> selected;

  // get input tracks
  const auto& input =
      ctx.eventStore.get<std::vector<VertexAndTracks>>(m_cfg.input);

  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  auto isValidTrack = [&](const auto& trk) {
    auto pos = trk.position(ctx.geoContext);
    auto dir = trk.unitDirection();
    auto rho = std::hypot(pos[Acts::eX], pos[Acts::eY]);
    auto phi = std::atan2(dir[Acts::eY], dir[Acts::eX]);
    auto eta = std::atanh(dir[Acts::eZ]);
    auto pt = trk.transverseMomentum();
    return within(rho, 0, m_cfg.rhoMax) and
           within(std::abs(pos[Acts::eZ]), 0, m_cfg.absZMax) and
           within(phi, m_cfg.phiMin, m_cfg.phiMax) and
           within(eta, m_cfg.etaMin, m_cfg.etaMax) and
           within(std::abs(eta), m_cfg.absEtaMin, m_cfg.absEtaMax) and
           within(pt, m_cfg.ptMin, m_cfg.ptMax) and
           (m_cfg.keepNeutral or (trk.charge() != 0));
  };

  for (const auto& vertexAndTracks : input) {
    VertexAndTracks sel;
    sel.vertex = vertexAndTracks.vertex;

    // Copy selected tracks over
    std::copy_if(vertexAndTracks.tracks.begin(), vertexAndTracks.tracks.end(),
                 std::back_inserter(sel.tracks), isValidTrack);
    // Only retain vertex if it still contains tracks
    if (not sel.tracks.empty()) {
      selected.push_back(std::move(sel));
    }
  }

  ACTS_DEBUG("event " << ctx.eventNumber << " selected " << selected.size()
                      << " from " << input.size() << " vertices.");

  // write selected tracks
  ctx.eventStore.add(m_cfg.output, std::move(selected));

  return ProcessCode::SUCCESS;
}
