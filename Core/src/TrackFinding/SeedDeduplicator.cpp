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

void SeedDeduplicator::reset(std::size_t nMeasKeys) {
  if (nMeasKeys >= kNone) {
    throw std::invalid_argument("SeedDeduplicator: too many keys");
  }
  m_trackHeads.assign(nMeasKeys, kNone);
  m_weights.assign(nMeasKeys, 1);
  m_nodes.clear();
  m_nTracks = 0;
}

SeedDeduplicator::Score SeedDeduplicator::duplicateThreshold(
    Score seedScore) const {
  const Score relative =
      seedScore > m_cfg.maxMissingScore ? seedScore - m_cfg.maxMissingScore : 0;
  return std::min(seedScore, std::max({Score{1}, m_cfg.minScore, relative}));
}

void SeedDeduplicator::addTrack(std::span<const MeasurementKey> measKeys) {
  const auto track = static_cast<std::uint32_t>(m_nTracks);
  for (const MeasurementKey measKey : measKeys) {
    assert(measKey < m_trackHeads.size() && "Key out of range");
    const std::uint32_t head = m_trackHeads[measKey];
    // a key that occurs twice on a track counts once
    if (head != kNone && m_nodes[head].track == track) {
      continue;
    }
    m_trackHeads[measKey] = static_cast<std::uint32_t>(m_nodes.size());
    m_nodes.push_back({track, head});
  }
  ++m_nTracks;
}

bool SeedDeduplicator::isDuplicate(
    std::span<const MeasurementKey> measKeys) const {
  Score seedScore = 0;
  for (const MeasurementKey measKey : measKeys) {
    seedScore += weight(measKey);
  }
  if (seedScore == 0) {
    return false;
  }
  const Score threshold = duplicateThreshold(seedScore);

  // the score that each track shares with the seed, for each track which holds
  // at least one seed key
  struct TrackScore {
    std::uint32_t track;
    Score sharedScore;
  };
  boost::container::small_vector<TrackScore, kInlineTracks> scores;
  Score bestSharedScore = 0;
  Score remaining = seedScore;
  for (const MeasurementKey measKey : measKeys) {
    const Weight w = weight(measKey);
    for (std::uint32_t node = m_trackHeads[measKey]; node != kNone;
         node = m_nodes[node].next) {
      const std::uint32_t track = m_nodes[node].track;
      auto it = std::ranges::find(scores, track, &TrackScore::track);
      if (it == scores.end()) {
        it = scores.insert(scores.end(), {track, 0});
      }
      it->sharedScore += w;
      if (it->sharedScore >= threshold) {
        return true;
      }
      bestSharedScore = std::max(bestSharedScore, it->sharedScore);
    }
    remaining -= w;
    // no track can reach the threshold with the keys that are left. This
    // check runs after each seed key, so it also stops the loop early for a
    // seed whose first keys no track holds.
    if (bestSharedScore + remaining < threshold) {
      return false;
    }
  }
  return false;
}

}  // namespace Acts
