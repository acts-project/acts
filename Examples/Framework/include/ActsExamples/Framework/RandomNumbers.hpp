// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/HashCombine.hpp"

#include <cstddef>
#include <cstdint>
#include <random>

#include <boost/functional/hash.hpp>

namespace ActsExamples {

struct AlgorithmContext;

/// The seed type used in the framework.
using RandomSeed = std::uint64_t;

/// The random number generator used in the framework.
///
/// Thin wrapper around std::mt19937_64 that remembers the seed it was
/// constructed with so that downstream consumers (e.g. lumi-block vertex
/// generators) can derive reproducible, user-seed-dependent sub-sequences.
class RandomEngine {
 public:
  using result_type = std::mt19937_64::result_type;

  /// Default constructor. Uses mt19937_64 default seed.
  RandomEngine() : m_seed(0) {}

  /// Construct with a specific seed.
  explicit RandomEngine(RandomSeed seed) : m_seed(seed), m_engine(seed) {}

  /// Return the seed this engine was constructed with.
  RandomSeed seed() const { return m_seed; }

  /// Create a new engine whose seed combines this engine's seed with an
  /// additional value. Useful for deriving deterministic sub-sequences
  /// (e.g. per-lumi-block seeds) that depend on the user's original seed.
  RandomEngine combinedWith(std::uint64_t extra) const {
    return RandomEngine(Acts::hashMixAndCombine(m_seed, extra));
  }

  result_type operator()() { return m_engine(); }

  static constexpr result_type min() { return std::mt19937_64::min(); }
  static constexpr result_type max() { return std::mt19937_64::max(); }

 private:
  RandomSeed m_seed;
  std::mt19937_64 m_engine;
};

/// Provide event and algorithm specific random number generator.s
///
/// This provides local random number generators, allowing for
/// thread-safe, lock-free, and reproducible random number generation across
/// single-threaded and multi-threaded test framework runs.
///
/// The role of the RandomNumbers is only to spawn local random number
/// generators. It does not, in and of itself, accommodate requests for specific
/// random number distributions (uniform, gaussian, etc). For this purpose,
/// clients should spawn their own local distribution objects
/// as needed, following the C++11 STL design.
class RandomNumbers {
 public:
  struct Config {
    RandomSeed seed = 1234567890u;  ///< random seed
  };

  explicit RandomNumbers(const Config& cfg);

  /// Spawn an algorithm-local random number generator. To avoid inefficiencies
  /// and multiple uses of a given RNG seed, this should only be done once per
  /// Algorithm invocation, after what the generator object should be reused.
  ///
  /// It calls generateSeed() for an event driven seed
  ///
  /// @param context is the AlgorithmContext of the host algorithm
  RandomEngine spawnGenerator(const AlgorithmContext& context) const;

  /// Generate a event and algorithm specific seed value.
  ///
  /// This should only be used in special cases e.g. where a custom
  /// random engine is used and `spawnGenerator` can not be used.
  RandomSeed generateSeed(const AlgorithmContext& context) const;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
