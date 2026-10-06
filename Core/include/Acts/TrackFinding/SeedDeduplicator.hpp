// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ranges>
#include <utility>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @brief Flags seeds whose measurements one accepted track already holds
///
/// Measurements are dense keys in `[0, nKeys)`. The caller chooses which keys
/// a seed contributes. Each key has a weight, 1 by default. A seed with the
/// summed weight `total` is a duplicate if one track holds keys of the summed
/// weight `min(total, max(1, minScore, total - maxMissingScore))`. The scores
/// of different tracks do not add up. Seeds can be queried in any order.
class SeedDeduplicator {
 public:
  /// Type of the dense measurement key
  using Key = std::uint32_t;
  /// Type of the weight of a key
  using Weight = std::uint8_t;
  /// Type of a summed weight
  using Score = std::uint32_t;

  /// Configuration of the duplicate threshold
  ///
  /// The default requires every key. A cut on `minScore` alone needs
  /// `maxMissingScore` at its maximum.
  struct Config {
    /// Minimum shared score that makes a seed a duplicate
    Score minScore = 0;
    /// Maximum score of the seed that the track may miss
    Score maxMissingScore = 0;
  };

  /// Constructor
  /// @param config The threshold configuration
  explicit SeedDeduplicator(const Config& config);

  /// Remove all tracks and set the key range. Every weight becomes 1.
  /// @param nKeys The number of keys, every key must be below this number
  void reset(std::size_t nKeys);

  /// @return The number of accepted tracks
  std::size_t nTracks() const { return m_nTracks; }

  /// Set the weight of a key
  /// @param key The key
  /// @param weight The weight of the key
  void setWeight(Key key, Weight weight) {
    assert(key < m_weights.size() && "Key out of range");
    m_weights[key] = weight;
  }

  /// @param key The key
  /// @return The weight of the key
  Weight weight(Key key) const {
    assert(key < m_weights.size() && "Key out of range");
    return m_weights[key];
  }

  /// @param total The summed weight of a seed
  /// @return The shared score which makes the seed a duplicate
  Score requiredScore(Score total) const;

  /// Accept a track
  /// @param keys The keys of the measurements on the track
  template <std::ranges::input_range range_t>
    requires std::convertible_to<std::ranges::range_value_t<range_t>, Key>
  void addTrack(range_t&& keys) {
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

  /// Check if one accepted track holds enough of the seed keys
  /// @param keys The keys of the seed, each key at most once
  /// @return True if the seed is a duplicate
  template <std::ranges::forward_range range_t>
    requires std::convertible_to<std::ranges::range_value_t<range_t>, Key>
  bool isDuplicate(range_t&& keys) const {
    Score remaining = 0;
    for (const Key key : keys) {
      remaining += weight(key);
    }
    if (remaining == 0) {
      return false;
    }
    const Score required = requiredScore(remaining);

    // the score of each track which holds at least one seed key
    boost::container::small_vector<std::pair<std::uint32_t, Score>, 8> scores;
    Score best = 0;
    for (const Key key : keys) {
      const Weight w = m_weights[key];
      for (std::uint32_t node = m_heads[key]; node != kNone;
           node = m_nodes[node].next) {
        const std::uint32_t track = m_nodes[node].track;
        auto it = std::ranges::find(scores, track,
                                    &std::pair<std::uint32_t, Score>::first);
        if (it == scores.end()) {
          it = scores.emplace(scores.end(), track, 0);
        }
        it->second += w;
        if (it->second >= required) {
          return true;
        }
        best = std::max(best, it->second);
      }
      remaining -= w;
      // no track can reach the threshold with the keys that are left
      if (best + remaining < required) {
        return false;
      }
    }
    return false;
  }

 private:
  static constexpr std::uint32_t kNone =
      std::numeric_limits<std::uint32_t>::max();

  /// One accepted track on one key
  struct Node {
    /// The index of the track
    std::uint32_t track;
    /// The next node of the same key, or `kNone`
    std::uint32_t next;
  };

  Config m_cfg;

  /// The latest node for each key, or `kNone`
  std::vector<std::uint32_t> m_heads;
  std::vector<Weight> m_weights;
  std::vector<Node> m_nodes;
  std::size_t m_nTracks = 0;
};

}  // namespace Acts
