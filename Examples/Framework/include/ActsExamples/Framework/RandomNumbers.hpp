// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <random>

#include <boost/container/small_vector.hpp>

namespace ActsExamples {
struct AlgorithmContext;

/// The random number generator used in the framework.
using RandomEngine = std::mt19937;  ///< Mersenne Twister

/// The seed type used in the framework.
using RandomSeed = std::uint32_t;

class RandomSeedSequence {
 public:
  using result_type = std::uint_least32_t;
  static constexpr std::size_t N = 4;
  using container_type = boost::container::small_vector<result_type, N>;

  RandomSeedSequence();
  explicit RandomSeedSequence(RandomSeed seed);
  template <typename InputIterator>
  RandomSeedSequence(InputIterator first, InputIterator last)
      : m_seed(first, last) {}
  template <typename InputRange>
  RandomSeedSequence(const InputRange& range)
      : m_seed(range.begin(), range.end()) {}
  RandomSeedSequence(std::initializer_list<result_type> il)
      : m_seed(il.begin(), il.end()) {}

  std::size_t size() const;

  container_type generate(std::size_t size) const;

  template <typename OutputIterator>
  void generate(OutputIterator first, OutputIterator last) const {
    container_type result = generate(static_cast<std::size_t>(last - first));
    std::copy(result.begin(), result.end(), first);
  }

  template <typename OutputIterator>
  void param(OutputIterator first) const {
    std::copy(m_seed.begin(), m_seed.end(), first);
  }

  RandomSeedSequence append(RandomSeed seed) const;

 private:
  container_type m_seed;
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

  RandomSeedSequence createEventAlgorithmSeedSequence(
      const AlgorithmContext& context) const;

  RandomEngine createEventAlgorithmEngine(
      const AlgorithmContext& context) const;

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
