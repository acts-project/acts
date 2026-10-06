// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/SeedDeduplicator.hpp"

#include <algorithm>
#include <stdexcept>

namespace Acts {

SeedDeduplicator::SeedDeduplicator(const Config& config) : m_cfg(config) {}

void SeedDeduplicator::reset(std::size_t nKeys) {
  if (nKeys >= kNone) {
    throw std::invalid_argument("SeedDeduplicator: too many keys");
  }
  m_heads.assign(nKeys, kNone);
  m_weights.assign(nKeys, 1);
  m_nodes.clear();
  m_nTracks = 0;
}

SeedDeduplicator::Score SeedDeduplicator::requiredScore(Score total) const {
  const Score relative =
      total > m_cfg.maxMissingScore ? total - m_cfg.maxMissingScore : 0;
  return std::min(total, std::max({Score{1}, m_cfg.minScore, relative}));
}

}  // namespace Acts
