// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Explicit random samplers, so that a seed reproduces an event on any
/// platform: `std::mt19937` is specified bit for bit, the standard
/// distributions on top of it are not.
///
/// @note Reproducible in the draws they take, not in the last bit of the
///       result: `std::log` and friends need not be exactly rounded.

#include <cmath>
#include <cstdint>
#include <numbers>
#include <optional>
#include <random>
#include <utility>

namespace ActsFatras::Synthetic {

struct SecondarySamplingConfig;

}  // namespace ActsFatras::Synthetic

namespace ActsFatras::Synthetic::detail {

/// Draw uniformly from [0, 1), with 24 bits of resolution.
inline float sampleUniform01(std::mt19937& rng) {
  // mt19937 yields 32 bits; keep the 24 a float can represent exactly
  return static_cast<float>(rng() >> 8) * 0x1.0p-24f;
}

/// Draw uniformly from (0, 1], i.e. what a logarithm can take: the grid of
/// `sampleUniform01` shifted up by one step.
inline float sampleUniform01Positive(std::mt19937& rng) {
  // the 24 bits run over [1, 2^24] instead of [0, 2^24), still exactly
  return static_cast<float>((rng() >> 8) + 1) * 0x1.0p-24f;
}

/// Draw a fair coin.
inline bool sampleBool(std::mt19937& rng) {
  // the top bit, which is the one `sampleUniform01(rng) < 0.5f` would look at
  return (rng() >> 31) != 0;
}

/// Draw either charge with equal probability, as plus or minus one.
inline float sampleCharge(std::mt19937& rng) {
  return sampleBool(rng) ? 1.f : -1.f;
}

/// Draw uniformly from `[lo, hi)`.
inline float sampleUniform(std::mt19937& rng, const float lo, const float hi) {
  return lo + (hi - lo) * sampleUniform01(rng);
}

/// Draw from a standard normal distribution, by Box-Muller. The second variate
/// is discarded rather than cached, so a call always costs two draws.
inline float sampleNormal(std::mt19937& rng) {
  constexpr float twoPi = 2.f * std::numbers::pi_v<float>;
  const float u1 = sampleUniform01Positive(rng);
  const float u2 = sampleUniform01(rng);
  return std::sqrt(-2.f * std::log(u1)) * std::cos(twoPi * u2);
}

/// Draw two independent standard normal variates, which Box-Muller makes at
/// the cost of one.
inline std::pair<float, float> sampleNormalPair(std::mt19937& rng) {
  constexpr float twoPi = 2.f * std::numbers::pi_v<float>;
  const float u1 = sampleUniform01Positive(rng);
  const float u2 = sampleUniform01(rng);
  const float radius = std::sqrt(-2.f * std::log(u1));
  return {radius * std::cos(twoPi * u2), radius * std::sin(twoPi * u2)};
}

/// Draw from a Rayleigh distribution, i.e. the length of a two-dimensional
/// Gaussian of width `scale` along either axis. That is neither the median nor
/// the mean of the result.
inline float sampleRayleigh(std::mt19937& rng, const float scale) {
  return scale * std::sqrt(-2.f * std::log(sampleUniform01Positive(rng)));
}

/// Draw a direction isotropically about an axis, as its cosine to that axis.
inline float sampleIsotropicCosine(std::mt19937& rng) {
  return sampleUniform(rng, -1.f, 1.f);
}

/// Draw from a log-normal distribution, `width` being the spread in e-folds.
inline float sampleLogNormal(std::mt19937& rng, const float median,
                             const float width) {
  return median * std::exp(width * sampleNormal(rng));
}

/// Draw from a Poisson distribution of non-negative `mean`, by Knuth's product
/// method: `mean + 1` draws on average, so fine for a small mean and slow for a
/// large one.
inline std::uint32_t samplePoisson(std::mt19937& rng, const float mean) {
  const float limit = std::exp(-mean);
  float product = sampleUniform01(rng);
  std::uint32_t count = 0;
  while (product > limit) {
    ++count;
    product *= sampleUniform01(rng);
  }
  return count;
}

/// Draw a transverse momentum from `dN/dpT ~ pT (1 + pT/ptScale)^-n` above
/// `minPt`, the Hagedorn spectrum in its invariant form, given a uniform `u` in
/// [0, 1). Its survival function does not invert, so five safeguarded Newton
/// steps solve it. The exponent `n` has to exceed two to normalise.
float samplePt(float u, float minPt, float ptScale, float n);

/// Draw a rapidity from a plateau with a Fermi edge, i.e. from
/// `dN/dy ~ 1 / (1 + exp((|y| - edge) / width))` over `[lo, hi)`, a `width` of
/// zero or less leaving it flat. By rejection, so only the shape moves and the
/// count the caller asked for is left alone.
float sampleRapidity(std::mt19937& rng, float lo, float hi, float edge,
                     float width);

/// The direction and momentum a secondary is produced with.
struct SecondaryKinematics {
  /// Azimuth of the momentum at the production point
  float phi{};
  /// Ratio of the longitudinal to the transverse momentum
  float cotTheta{};
  /// Transverse momentum
  float pt{};
};

/// Draw where a secondary goes and how hard it is, the channel deciding both
/// the laws and the axis it comes off: the parent's direction (`trackPhi`,
/// `trackCotTheta`) or the line from the beam line through the production point
/// (`radialPhi`, `radialCotTheta`) that an electron follows. `electronShare` is
/// the chance of that electron channel, `parentP` the momentum the daughter is
/// bounded by. Nothing below `SecondarySamplingConfig::minPt`.
std::optional<SecondaryKinematics> sampleSecondary(
    std::mt19937& rng, const SecondarySamplingConfig& cfg, float trackPhi,
    float trackCotTheta, float radialPhi, float radialCotTheta, float parentP,
    float electronShare);

}  // namespace ActsFatras::Synthetic::detail
