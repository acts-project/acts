// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/TruthTracking/TrackSelector.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/TruthTracking/TruthVerticesToTracks.hpp"
#include "Acts/EventData/TrackParameters.hpp"

FW::TrackSelector::TrackSelector(const Config& cfg, Acts::Logging::Level level)
    : FW::BareAlgorithm("Selector", level), m_cfg(cfg) {
  if (m_cfg.input.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
  if (m_cfg.output.empty()) {
    throw std::invalid_argument("Missing output collection");
  }
}

FW::ProcessCode FW::TrackSelector::execute(
    const FW::AlgorithmContext& ctx) const {
  std::vector<VertexAndTracks> selected;

  // get input tracks
  const auto& input =
      ctx.eventStore.get<std::vector<VertexAndTracks>>(m_cfg.input);

  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  auto isValidTrack = [&](const auto& trk) {
    auto rho = std::hypot(trk.position().x(), trk.position().y());
    auto phi = std::atan2(trk.momentum().y(), trk.momentum().x());
    auto eta = std::atanh(trk.momentum().z() / trk.momentum().norm());
    auto pt = std::hypot(trk.momentum().x(), trk.momentum().y());
    return within(rho, 0, m_cfg.rhoMax) and
           within(std::abs(trk.position().z()), 0, m_cfg.absZMax) and
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
