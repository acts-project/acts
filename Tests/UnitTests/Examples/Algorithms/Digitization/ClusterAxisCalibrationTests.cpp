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
#include "ActsExamples/Digitization/ClusterAxisCalibration.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <random>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

namespace ActsTests {

namespace {

/// A sensor with equidistant bins on one axis, from -edge to +edge
struct AxisSensor {
  double pitch = 50_um;
  double thickness = 150_um;
  int halfBins = 50;
  ClusterAxisReadout readout;

  /// The footprint at the angle `alpha` to the normal. The charge is the path
  /// length, as in the geometric digitization.
  ClusterAxisFootprint footprint(double alpha) const {
    return {thickness * std::tan(alpha), thickness / std::cos(alpha)};
  }

  /// Deposit a track with the midpoint `mid` and add the cell noise. The
  /// result maps the fired bins to their charge.
  std::map<int, double> deposit(double mid, const ClusterAxisFootprint& f,
                                std::mt19937& rng) const {
    std::normal_distribution<double> gauss(0., 1.);
    std::map<int, double> fired;
    const double low = mid - 0.5 * f.length;
    const double high = mid + 0.5 * f.length;
    const int first =
        std::max(static_cast<int>(std::floor(low / pitch)), -halfBins);
    const int last =
        std::min(static_cast<int>(std::floor(high / pitch)), halfBins - 1);
    for (int bin = first; bin <= last; ++bin) {
      const double overlap =
          std::min(high, (bin + 1) * pitch) - std::max(low, bin * pitch);
      double charge =
          f.length > 0 ? f.pathLength * overlap / f.length : f.pathLength;
      if (charge <= 0) {
        continue;
      }
      charge = std::max(0., charge + readout.noise * gauss(rng));
      if (charge > readout.threshold) {
        fired[bin] = charge;
      }
    }
    return fired;
  }

  ClusterAxisProjection project(const std::map<int, double>& fired) const {
    const int first = fired.begin()->first;
    const int last = fired.rbegin()->first;
    ClusterAxisProjection p;
    p.firstLower = first * pitch;
    p.firstUpper = (first + 1) * pitch;
    p.lastLower = last * pitch;
    p.lastUpper = (last + 1) * pitch;
    p.firstCharge = fired.begin()->second;
    p.lastCharge = fired.rbegin()->second;
    p.beforeWidth = first > -halfBins ? pitch : 0;
    p.afterWidth = last < halfBins - 1 ? pitch : 0;
    return p;
  }
};

}  // namespace

BOOST_AUTO_TEST_SUITE(DigitizationSuite)

// The calibration works in any charge unit, so a charge calibration can run
// before it, for example from the time over threshold.
BOOST_AUTO_TEST_CASE(ClusterAxisChargeUnit) {
  AxisSensor sensor;
  sensor.readout = {5_um, 5_um, true};
  const auto f = sensor.footprint(30_degree);
  std::mt19937 rng(7);
  const auto p = sensor.project(sensor.deposit(0.3 * sensor.pitch, f, rng));
  const auto reference = calibrateClusterAxis(p, f, sensor.readout);

  // About 80 electrons per micrometre
  const double scale = 8e4;
  auto scaled = p;
  scaled.firstCharge *= scale;
  scaled.lastCharge *= scale;
  const ClusterAxisReadout scaledReadout{sensor.readout.threshold * scale,
                                         sensor.readout.noise * scale, true,
                                         scale};
  const auto result = calibrateClusterAxis(scaled, f, scaledReadout);
  CHECK_CLOSE_REL(result.mean, reference.mean, 1e-9);
  CHECK_CLOSE_REL(result.variance, reference.variance, 1e-9);
}

// The posterior describes the clusters of the model: the pulls have unit
// width for every readout
BOOST_DATA_TEST_CASE(ClusterAxisPulls,
                     boost::unit_test::data::make({5., 20., 40.}) *
                         boost::unit_test::data::make({false, true}) *
                         boost::unit_test::data::make({2., 10.}),
                     angle, analogue, noise) {
  AxisSensor sensor;
  const auto f = sensor.footprint(angle * 1_degree);
  sensor.readout = {0.05 * f.pathLength, noise * 1_um, analogue};
  std::mt19937 rng(1234);
  std::uniform_real_distribution<double> uniform(0., sensor.pitch);

  double sumPull2 = 0;
  std::size_t count = 0;
  for (std::size_t i = 0; i < 2000; ++i) {
    const double mid = uniform(rng);
    const auto fired = sensor.deposit(mid, f, rng);
    if (fired.empty()) {
      continue;
    }
    const auto result =
        calibrateClusterAxis(sensor.project(fired), f, sensor.readout);
    const double pull = (result.mean - mid) / std::sqrt(result.variance);
    if (std::abs(pull) < 5) {
      sumPull2 += pull * pull;
      ++count;
    }
  }
  BOOST_REQUIRE_GT(count, 1900u);
  const double pullRms = std::sqrt(sumPull2 / static_cast<double>(count));
  BOOST_TEST_MESSAGE("angle " << angle << " analogue " << analogue << " noise "
                              << noise << ": pull rms " << pullRms);
  CHECK_CLOSE_REL(pullRms, 1., 0.05);
}

// A cluster that the footprint cannot make is uniform over the cluster
BOOST_AUTO_TEST_CASE(ClusterAxisFallback) {
  AxisSensor sensor;
  const auto f = sensor.footprint(10_degree);
  ClusterAxisProjection p;
  p.firstLower = 0;
  p.firstUpper = sensor.pitch;
  p.lastLower = 3 * sensor.pitch;
  p.lastUpper = 4 * sensor.pitch;
  p.beforeWidth = p.afterWidth = sensor.pitch;
  const auto result = calibrateClusterAxis(p, f, sensor.readout);
  CHECK_CLOSE_ABS(result.mean, 2 * sensor.pitch, 1e-12);
  CHECK_CLOSE_REL(result.variance, std::pow(4 * sensor.pitch, 2) / 12, 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
