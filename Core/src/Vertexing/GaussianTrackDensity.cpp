// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/GaussianTrackDensity.hpp"

double Acts::GaussianTrackDensity::globalMaximum(
    const std::vector<const Acts::BoundParameters>& trackList,
    State& state) const {
  addTracks(trackList, state);
  return state.trackDensity.globalMaximum(state.trackDensityState);
}

std::pair<double, double> Acts::GaussianTrackDensity::globalMaximumWithWidth(
    const std::vector<const Acts::BoundParameters>& trackList,
    State& state) const {
  addTracks(trackList, state);
  return state.trackDensity.globalMaximumWithWidth(state.trackDensityState);
}

void Acts::GaussianTrackDensity::addTracks(
    const std::vector<const Acts::BoundParameters>& trackList,
    State& state) const {
  const double d0SignificanceCut =
      m_cfg.d0MaxSignificance * m_cfg.d0MaxSignificance;
  const double z0SignificanceCut =
      m_cfg.z0MaxSignificance * m_cfg.z0MaxSignificance;

  for (const auto& trk : trackList) {
    state.trackDensity.addTrack(state.trackDensityState, trk, d0SignificanceCut,
                                z0SignificanceCut);
  }
}