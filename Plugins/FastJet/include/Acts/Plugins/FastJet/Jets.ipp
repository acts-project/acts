// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

namespace Acts::FastJet {

template <typename TrackContainer>
std::vector<fastjet::PseudoJet> InputTracks<TrackContainer>::fourMomenta()
    const {
  std::vector<fastjet::PseudoJet> inputs;
  for (std::size_t i = 0; i < m_tracks.size(); i++) {
    Acts::Vector4 p = m_tracks.getTrack(i).fourMomentum();
    inputs.emplace_back(p[Acts::eMom0], p[Acts::eMom1], p[Acts::eMom2],
                        p[Acts::eEnergy]);
    inputs.back().set_user_index(i);
  }
  return inputs;
}

template <typename TrackContainer>
std::vector<typename TrackContainer::TrackProxy>
InputTracks<TrackContainer>::tracksInJet(const fastjet::PseudoJet& jet,
                                         std::optional<float> coreR) {
  fastjet::Selector sel = fastjet::SelectorIdentity();
  if (coreR.has_value()) {
    if (*coreR < 0) {
      throw std::invalid_argument("coreR must be positive!");
    }
    sel = fastjet::SelectorCircle(*coreR);
    sel.set_reference(jet);
  }

  std::vector<typename TrackContainer::TrackProxy> tracks;
  for (fastjet::PseudoJet& cst : sel(jet.constituents())) {
    tracks.push_back(m_tracks.getTrack(cst.user_index()));
  }

  return tracks;
}

}  // namespace Acts::FastJet
