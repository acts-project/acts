// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Digitization/ClusterOnTrackCalibrator.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <random>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;
using namespace ActsExamples::detail;

namespace ActsTests {

namespace {

/// A footprint at the angle `alpha` to the normal of a sensor of thickness
/// `thickness`, without a direction error.
ClusterFootprint makeFootprint(double thickness, double alpha) {
  ClusterFootprint footprint;
  footprint.length = thickness * std::abs(std::tan(alpha));
  footprint.pathLength = thickness / std::cos(alpha);
  return footprint;
}

/// Deposit the footprint with the midpoint `mid` on bins of width `pitch`,
/// the same way the digitization does, and project the cluster.
std::optional<ClusterProjection> makeCluster(double mid, double pitch,
                                             const ClusterFootprint& footprint,
                                             double threshold,
                                             double chargeNoise,
                                             std::mt19937& rng) {
  const double low = mid - 0.5 * footprint.length;
  const double high = mid + 0.5 * footprint.length;
  std::normal_distribution<double> noise(0., chargeNoise);
  std::map<int, double> fired;
  for (int bin = static_cast<int>(std::floor(low / pitch));
       bin <= static_cast<int>(std::floor(high / pitch)); ++bin) {
    const double overlap =
        std::min(high, (bin + 1) * pitch) - std::max(low, bin * pitch);
    double charge = footprint.length > 0
                        ? overlap / footprint.length * footprint.pathLength
                        : footprint.pathLength;
    if (chargeNoise > 0) {
      charge = std::max(0., charge + noise(rng));
    }
    if (charge > threshold) {
      fired[bin] = charge;
    }
  }
  if (fired.empty()) {
    return std::nullopt;
  }
  const auto& [first, firstCharge] = *fired.begin();
  const auto& [last, lastCharge] = *fired.rbegin();
  ClusterProjection projection;
  projection.firstLower = first * pitch;
  projection.firstUpper = (first + 1) * pitch;
  projection.lastLower = last * pitch;
  projection.lastUpper = (last + 1) * pitch;
  projection.firstCharge = firstCharge;
  projection.lastCharge = lastCharge;
  projection.firstCells = 1;
  projection.lastCells = 1;
  projection.beforeWidth = pitch;
  projection.afterWidth = pitch;
  return projection;
}

ClusterProjection makeBins(double pitch, int first, int last) {
  ClusterProjection projection;
  projection.firstLower = first * pitch;
  projection.firstUpper = (first + 1) * pitch;
  projection.lastLower = last * pitch;
  projection.lastUpper = (last + 1) * pitch;
  projection.firstCells = 1;
  projection.lastCells = 1;
  return projection;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(DigitizationSuite)

// Pitch 50 um, thickness 150 um, 10 degrees: the footprint is 0.529 pitches.
BOOST_AUTO_TEST_CASE(DigitalBoxWidths) {
  const double pitch = 50_um;
  const ClusterFootprint footprint = makeFootprint(150_um, 10_degree);
  const double phi = footprint.length / pitch;

  const auto single =
      calibrateClusterAxis(makeBins(pitch, 2, 2), footprint, 0, false, 0);
  BOOST_CHECK(single.consistent);
  CHECK_CLOSE_ABS(single.position, 2.5 * pitch, 1e-9);
  CHECK_CLOSE_REL(std::sqrt(single.variance), (1 - phi) * pitch / std::sqrt(12),
                  1e-9);

  const auto pair =
      calibrateClusterAxis(makeBins(pitch, 2, 3), footprint, 0, false, 0);
  BOOST_CHECK(pair.consistent);
  CHECK_CLOSE_ABS(pair.position, 3 * pitch, 1e-9);
  CHECK_CLOSE_REL(std::sqrt(pair.variance), phi * pitch / std::sqrt(12), 1e-9);

  // Three bins cannot come from 0.529 pitches: uniform over the cluster
  const auto triple =
      calibrateClusterAxis(makeBins(pitch, 2, 4), footprint, 0, false, 0);
  BOOST_CHECK(!triple.consistent);
  CHECK_CLOSE_ABS(triple.position, 3.5 * pitch, 1e-9);
  CHECK_CLOSE_REL(triple.variance, std::pow(3 * pitch, 2) / 12, 1e-9);
}

// The threshold shrinks the footprint by the threshold fraction at both ends
BOOST_AUTO_TEST_CASE(DigitalThreshold) {
  const double pitch = 50_um;
  const ClusterFootprint footprint = makeFootprint(150_um, 10_degree);
  const double threshold = 0.08 * footprint.pathLength;
  const double phi = footprint.length * (1 - 2 * 0.08) / pitch;

  const auto single = calibrateClusterAxis(makeBins(pitch, 0, 0), footprint,
                                           threshold, false, 0);
  CHECK_CLOSE_REL(std::sqrt(single.variance), (1 - phi) * pitch / std::sqrt(12),
                  1e-9);
  const auto pair = calibrateClusterAxis(makeBins(pitch, 0, 1), footprint,
                                         threshold, false, 0);
  CHECK_CLOSE_REL(std::sqrt(pair.variance), phi * pitch / std::sqrt(12), 1e-9);
}

// Without noise, the charge of the edge bins gives the exact position
BOOST_AUTO_TEST_CASE(AnalogueExact) {
  const double pitch = 50_um;
  const ClusterFootprint footprint = makeFootprint(150_um, 25_degree);
  std::mt19937 rng(42);
  const double mid = 3.3 * pitch;
  const auto cluster = makeCluster(mid, pitch, footprint, 0, 0, rng);
  BOOST_REQUIRE(cluster.has_value());
  const auto result = calibrateClusterAxis(*cluster, footprint, 0, true, 0);
  BOOST_CHECK(result.consistent);
  CHECK_CLOSE_ABS(result.position, mid, 1e-9);
  CHECK_SMALL(result.variance, 1e-12);
}

// Without noise, the truth is always in the digital box. The pulls have unit
// width for every readout.
BOOST_DATA_TEST_CASE(Pulls,
                     boost::unit_test::data::make({0., 5., 10., 20., 40.}) *
                         boost::unit_test::data::make({false, true}) *
                         boost::unit_test::data::make({0., 5., 10.}),
                     angle, analogue, noise) {
  const double pitch = 50_um;
  const double thickness = 150_um;
  const double chargeNoise = noise * 1_um;
  const ClusterFootprint footprint = makeFootprint(thickness, angle * 1_degree);
  const double threshold = 0.05 * footprint.pathLength;

  std::mt19937 rng(1234);
  std::uniform_real_distribution<double> uniform(0., pitch);

  double sumPull2 = 0;
  // The clusters that the footprint cannot make
  std::size_t n = 0;
  std::size_t consistent = 0;
  for (std::size_t i = 0; i < 20000; ++i) {
    const double mid = 10 * pitch + uniform(rng);
    const auto cluster =
        makeCluster(mid, pitch, footprint, threshold, chargeNoise, rng);
    if (!cluster.has_value()) {
      continue;
    }
    const auto result = calibrateClusterAxis(*cluster, footprint, threshold,
                                             analogue, chargeNoise);
    if (!result.consistent) {
      ++n;
      continue;
    }
    if (result.variance == 0) {
      // Analogue readout without noise is exact, unless the threshold
      // removed the cells after both edge bins
      CHECK_CLOSE_ABS(result.position, mid, 1e-9);
      continue;
    }
    ++consistent;
    if (!analogue && chargeNoise == 0) {
      const double half = std::sqrt(3 * result.variance);
      BOOST_CHECK_LE(std::abs(mid - result.position), half + 1e-12);
    }
    sumPull2 += std::pow(mid - result.position, 2) / result.variance;
  }
  BOOST_TEST_MESSAGE("angle " << angle << " analogue " << analogue << " noise "
                              << noise << " consistent " << consistent
                              << " inconsistent " << n << " pull rms "
                              << std::sqrt(sumPull2 / consistent));
  BOOST_CHECK_LE(n, 0.01 * consistent);
  if (consistent > 0) {
    CHECK_CLOSE_REL(std::sqrt(sumPull2 / consistent), 1., 0.05);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
