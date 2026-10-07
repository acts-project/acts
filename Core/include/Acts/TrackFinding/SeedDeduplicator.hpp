// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

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
  /// @return The score that one track must share with the seed to make the
  ///         seed a duplicate
  Score duplicateThreshold(Score total) const;

  /// Accept a track
  /// @param keys The keys of the measurements on the track
  void addTrack(std::span<const Key> keys);

  /// Check if one accepted track holds enough of the seed keys
  /// @param keys The keys of the seed, each key at most once
  /// @return True if the seed is a duplicate
  bool isDuplicate(std::span<const Key> keys) const;

 private:
  static constexpr std::uint32_t kNone =
      std::numeric_limits<std::uint32_t>::max();
  /// Number of tracks per seed query that need no heap allocation
  static constexpr std::size_t kInlineTracks = 8;

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
