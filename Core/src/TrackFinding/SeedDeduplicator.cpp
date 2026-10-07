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

#include <boost/container/small_vector.hpp>

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

SeedDeduplicator::Score SeedDeduplicator::duplicateThreshold(
    Score total) const {
  const Score relative =
      total > m_cfg.maxMissingScore ? total - m_cfg.maxMissingScore : 0;
  return std::min(total, std::max({Score{1}, m_cfg.minScore, relative}));
}

void SeedDeduplicator::addTrack(std::span<const Key> keys) {
  const auto track = static_cast<std::uint32_t>(m_nTracks);
  for (const Key key : keys) {
    assert(key < m_heads.size() && "Key out of range");
    const std::uint32_t head = m_heads[key];
    // a key that occurs twice on a track counts once
    if (head != kNone && m_nodes[head].track == track) {
      continue;
    }
    m_heads[key] = static_cast<std::uint32_t>(m_nodes.size());
    m_nodes.push_back({track, head});
  }
  ++m_nTracks;
}

bool SeedDeduplicator::isDuplicate(std::span<const Key> keys) const {
  Score remaining = 0;
  for (const Key key : keys) {
    remaining += weight(key);
  }
  if (remaining == 0) {
    return false;
  }
  const Score threshold = duplicateThreshold(remaining);

  // the score of each track which holds at least one seed key
  struct TrackScore {
    std::uint32_t track;
    Score score;
  };
  boost::container::small_vector<TrackScore, kInlineTracks> scores;
  Score best = 0;
  for (const Key key : keys) {
    const Weight w = m_weights[key];
    for (std::uint32_t node = m_heads[key]; node != kNone;
         node = m_nodes[node].next) {
      const std::uint32_t track = m_nodes[node].track;
      auto it = std::ranges::find(scores, track, &TrackScore::track);
      if (it == scores.end()) {
        it = scores.insert(scores.end(), {track, 0});
      }
      it->score += w;
      if (it->score >= threshold) {
        return true;
      }
      best = std::max(best, it->score);
    }
    remaining -= w;
    // no track can reach the threshold with the keys that are left
    if (best + remaining < threshold) {
      return false;
    }
  }
  return false;
}

}  // namespace Acts
