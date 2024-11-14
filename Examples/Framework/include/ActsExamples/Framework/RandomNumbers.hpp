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
#include <type_traits>

namespace ActsExamples {
struct AlgorithmContext;

/// The random number generator used in the framework.
using RandomEngine = std::mt19937;  ///< Mersenne Twister

/// The seed type used in the framework.
using RandomSeed = std::uint32_t;

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

/// Distribute entropy from a sequence of input values to a sequence of output
/// values.
///
/// This function is copied from std::seed_seq::generate() and modified to avoid
/// the intermediate std::vector.
///
/// @tparam InputIteratorBegin is the type of the input iterator begin
/// @tparam InputIteratorEnd is the type of the input iterator end
/// @tparam OutputIteratorBegin is the type of the output iterator begin
/// @tparam OutputIteratorEnd is the type of the output iterator end
///
/// @param inputIteratorBegin is the begin of the input sequence
/// @param inputIteratorEnd is the end of the input sequence
/// @param outputIteratorBegin is the begin of the output sequence
/// @param outputIteratorEnd is the end of the output sequence
template <typename InputIteratorBegin, typename InputIteratorEnd,
          typename OutputIteratorBegin, typename OutputIteratorEnd>
void distributeEntropy(InputIteratorBegin inputIteratorBegin,
                       InputIteratorEnd inputIteratorEnd,
                       OutputIteratorBegin outputIteratorBegin,
                       OutputIteratorEnd outputIteratorEnd)
  requires((std::input_iterator<InputIteratorBegin> &&
            std::random_access_iterator<InputIteratorBegin> &&
            std::is_assignable_v<std::uint_least32_t,
                                 decltype(*inputIteratorBegin)>) &&
           (std::input_iterator<InputIteratorEnd> &&
            std::random_access_iterator<InputIteratorEnd> &&
            std::is_assignable_v<std::uint_least32_t,
                                 decltype(*inputIteratorEnd)>) &&
           (std::output_iterator<OutputIteratorBegin, std::uint_least32_t> &&
            std::random_access_iterator<OutputIteratorBegin>) &&
           (std::output_iterator<OutputIteratorEnd, std::uint_least32_t> &&
            std::random_access_iterator<OutputIteratorEnd>))
{
  using result_type = std::uint_least32_t;

  if (inputIteratorBegin == inputIteratorEnd) {
    return;
  }
  if (outputIteratorBegin == outputIteratorEnd) {
    throw std::invalid_argument("outputIteratorBegin == outputIteratorEnd");
  }

  const auto __v_ = inputIteratorBegin;
  const auto __v_size =
      static_cast<std::size_t>(inputIteratorEnd - inputIteratorBegin);
  const auto __first = outputIteratorBegin;
  const auto __last = outputIteratorEnd;

  std::fill(__first, __last, 0x8b8b8b8b);
  const std::size_t __n = static_cast<std::size_t>(__last - __first);
  const std::size_t __s = __v_size;
  const std::size_t __t = (__n >= 623)  ? 11
                          : (__n >= 68) ? 7
                          : (__n >= 39) ? 5
                          : (__n >= 7)  ? 3
                                        : (__n - 1) / 2;
  const std::size_t __p = (__n - __t) / 2;
  const std::size_t __q = __p + __t;
  const std::size_t __m = std::max(__s + 1, __n);
  // __k = 0;
  {
    result_type __r =
        1664525 * _Tp(__first[0] ^ __first[__p] ^ __first[__n - 1]);
    __first[__p] += __r;
    __r += __s;
    __first[__q] += __r;
    __first[0] = __r;
  }
  // Initialize indexing terms used with if statements as an optimization to
  // avoid calculating modulo n on every loop iteration for each term.
  std::size_t __kmodn = 0;           // __k % __n
  std::size_t __k1modn = __n - 1;    // (__k - 1) % __n
  std::size_t __kpmodn = __p % __n;  // (__k + __p) % __n
  std::size_t __kqmodn = __q % __n;  // (__k + __q) % __n

  for (std::size_t __k = 1; __k <= __s; ++__k) {
    if (++__kmodn == __n) {
      __kmodn = 0;
    }
    if (++__k1modn == __n) {
      __k1modn = 0;
    }
    if (++__kpmodn == __n) {
      __kpmodn = 0;
    }
    if (++__kqmodn == __n) {
      __kqmodn = 0;
    }

    result_type __r =
        1664525 * _Tp(__first[__kmodn] ^ __first[__kpmodn] ^ __first[__k1modn]);
    __first[__kpmodn] += __r;
    __r += __kmodn + __v_[__k - 1];
    __first[__kqmodn] += __r;
    __first[__kmodn] = __r;
  }
  for (std::size_t __k = __s + 1; __k < __m; ++__k) {
    if (++__kmodn == __n) {
      __kmodn = 0;
    }
    if (++__k1modn == __n) {
      __k1modn = 0;
    }
    if (++__kpmodn == __n) {
      __kpmodn = 0;
    }
    if (++__kqmodn == __n) {
      __kqmodn = 0;
    }

    result_type __r =
        1664525 * _Tp(__first[__kmodn] ^ __first[__kpmodn] ^ __first[__k1modn]);
    __first[__kpmodn] += __r;
    __r += __kmodn;
    __first[__kqmodn] += __r;
    __first[__kmodn] = __r;
  }
  for (std::size_t __k = __m; __k < __m + __n; ++__k) {
    if (++__kmodn == __n) {
      __kmodn = 0;
    }
    if (++__k1modn == __n) {
      __k1modn = 0;
    }
    if (++__kpmodn == __n) {
      __kpmodn = 0;
    }
    if (++__kqmodn == __n) {
      __kqmodn = 0;
    }

    result_type __r = 1566083941 * _Tp(__first[__kmodn] + __first[__kpmodn] +
                                       __first[__k1modn]);
    __first[__kpmodn] ^= __r;
    __r -= __kmodn;
    __first[__kqmodn] ^= __r;
    __first[__kmodn] = __r;
  }
}

}  // namespace ActsExamples
