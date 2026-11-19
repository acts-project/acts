// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/ClusterAxisCalibration.hpp"

#include "Acts/Utilities/AlgebraHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>

namespace ActsExamples {

namespace {

constexpr double kMinusInfinity = -std::numeric_limits<double>::infinity();

/// The best log likelihood below which the cluster counts as incompatible
/// with the footprint. This is a chi2 of 25.
constexpr double kMinLogLikelihood = -12.5;

/// The log of the standard normal CDF, without underflow in the tails. Both
/// branches get clamped arguments, because the compiler can evaluate both.
double logNormalCdf(double x) {
  constexpr double kTail = 8.;
  const double inner = std::clamp(x, -kTail, kTail);
  const double tail = std::min(x, -kTail);
  const double tail2 = tail * tail;
  // Asymptotic expansion of the Mills ratio
  const double asymptotic = -0.5 * tail2 - std::log(-tail) -
                            0.5 * std::log(2. * std::numbers::pi) +
                            std::log1p(-1. / tail2 + 3. / (tail2 * tail2));
  const double central =
      std::log(0.5 * std::erfc(-inner / std::numbers::sqrt2));
  return x > kTail ? 0. : (x < -kTail ? asymptotic : central);
}

/// The fraction of the segment with midpoint `mid` and length `length` in
/// the bin [lower, upper). The division is safe for a zero length, because
/// the compiler can evaluate it before the check.
double overlapFraction(double mid, double length, double lower, double upper) {
  const double overlap =
      std::min(upper, mid + 0.5 * length) - std::max(lower, mid - 0.5 * length);
  const double fraction = std::max(overlap, 0.) /
                          std::max(length, std::numeric_limits<double>::min());
  const double inside = (lower <= mid && mid < upper) ? 1. : 0.;
  return length > 0 ? fraction : inside;
}

}  // namespace

double clusterAxisLogLikelihood(const ClusterAxisProjection& projection,
                                const ClusterAxisFootprint& footprint,
                                const ClusterAxisReadout& readout, double mid) {
  const auto& p = projection;
  const double length = footprint.length;
  const double pathCharge = readout.chargePerLength * footprint.pathLength;
  // Without noise, a small noise keeps the likelihood finite
  const double noise = std::max(
      {readout.noise, 1e-9 * pathCharge, std::numeric_limits<double>::min()});
  const auto charge = [&](double lower, double upper) {
    return pathCharge * overlapFraction(mid, length, lower, upper);
  };
  const auto logFire = [&](double deposit, double sigma, bool fires) {
    const double z = (deposit - readout.threshold) / sigma;
    return logNormalCdf(fires ? z : -z);
  };
  const auto logFired = [&](double deposit, double observed,
                            std::size_t cells) {
    const double sigma = noise * std::sqrt(static_cast<double>(cells));
    if (!readout.analogue) {
      return logFire(deposit, sigma, true);
    }
    const double z = (observed - deposit) / sigma;
    return -0.5 * z * z;
  };

  const double first = charge(p.firstLower, p.firstUpper);
  const double last = charge(p.lastLower, p.lastUpper);
  // The segment touches every fired bin
  if (first <= 0 || last <= 0) {
    return kMinusInfinity;
  }
  double result = logFired(first, p.firstCharge, p.firstCells);
  if (p.firstLower != p.lastLower) {
    result += logFired(last, p.lastCharge, p.lastCells);
  }
  // A touched neighbour stays below the threshold, an untouched one never
  // fires
  if (p.beforeWidth > 0) {
    const double before = charge(p.firstLower - p.beforeWidth, p.firstLower);
    if (before > 0) {
      result += logFire(before, noise, false);
    }
  }
  if (p.afterWidth > 0) {
    const double after = charge(p.lastUpper, p.lastUpper + p.afterWidth);
    if (after > 0) {
      result += logFire(after, noise, false);
    }
  }
  return result;
}

ClusterAxisPosterior calibrateClusterAxis(
    const ClusterAxisProjection& projection,
    const ClusterAxisFootprint& footprint, const ClusterAxisReadout& readout) {
  constexpr std::size_t kPoints = 64;
  constexpr std::size_t kMaxPasses = 8;
  /// The log likelihood below the maximum that counts as outside the support
  constexpr double kSupport = 20;

  const auto& p = projection;
  // Uniform over the whole cluster if the footprint cannot make the cluster
  const ClusterAxisPosterior fallback{
      0.5 * (p.firstLower + p.lastUpper),
      std::pow(p.lastUpper - p.firstLower, 2) / 12.};

  // The positions where the segment touches both edge bins. The midpoint is
  // on the sensor, even when the footprint leaves it.
  double low = p.lastLower - 0.5 * footprint.length;
  double high = p.firstUpper + 0.5 * footprint.length;
  if (p.beforeWidth == 0) {
    low = std::max(low, p.firstLower);
  }
  if (p.afterWidth == 0) {
    high = std::min(high, p.lastUpper);
  }
  if (!(footprint.pathLength > 0) || !(low < high)) {
    return fallback;
  }

  std::array<double, kPoints> logL{};
  ClusterAxisPosterior result;
  // The grid can miss a narrow peak, so the best value of all passes counts
  double bestLog = kMinusInfinity;
  for (std::size_t pass = 0; pass < kMaxPasses; ++pass) {
    const double step = (high - low) / kPoints;
    const auto position = [&](std::size_t i) {
      return low + (static_cast<double>(i) + 0.5) * step;
    };
    double maxLog = kMinusInfinity;
    for (std::size_t i = 0; i < kPoints; ++i) {
      logL[i] = clusterAxisLogLikelihood(p, footprint, readout, position(i));
      maxLog = std::max(maxLog, logL[i]);
    }
    if (maxLog == kMinusInfinity) {
      return fallback;
    }
    bestLog = std::max(bestLog, maxLog);

    double norm = 0;
    double sum1 = 0;
    double sum2 = 0;
    std::size_t first = kPoints;
    std::size_t last = 0;
    const double origin = 0.5 * (low + high);
    for (std::size_t i = 0; i < kPoints; ++i) {
      const double w = Acts::safeExp(logL[i] - maxLog);
      const double delta = position(i) - origin;
      norm += w;
      sum1 += w * delta;
      sum2 += w * delta * delta;
      if (logL[i] > maxLog - kSupport) {
        first = std::min(first, i);
        last = i;
      }
    }
    const double shift = sum1 / norm;
    result.mean = origin + shift;
    // The grid spreads the weight of a point uniformly over its step
    result.variance =
        std::max(sum2 / norm - shift * shift, 0.) + step * step / 12.;

    // Refine the grid to the support, if that is a part of the grid
    const double supportLow = position(first) - step;
    const double supportHigh = position(last) + step;
    if (supportHigh - supportLow > 0.5 * (high - low)) {
      break;
    }
    low = std::max(low, supportLow);
    high = std::min(high, supportHigh);
  }
  return bestLog > kMinLogLikelihood ? result : fallback;
}

}  // namespace ActsExamples
