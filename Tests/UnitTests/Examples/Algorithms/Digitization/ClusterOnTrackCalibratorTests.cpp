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
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"
#include "ActsExamples/Digitization/ClusterOnTrackCalibrator.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsFatras/Digitization/Channelizer.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <random>
#include <utility>

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace ActsExamples;

namespace ActsTests {

namespace {

constexpr int kBinOffset = 50;

/// A sensor with one measured axis. Its normal is the global x axis, so a
/// track at the angle `alpha` to the normal in the local x-z plane has
/// phi = alpha and theta = 90 degrees.
struct Sensor {
  double pitch = 50_um;
  double thickness = 150_um;
  double threshold = 0;
  bool analogue = false;
  double chargeNoise = 0;

  GeometryIdentifier geoId =
      GeometryIdentifier().withVolume(1).withLayer(2).withSensitive(3);
  std::shared_ptr<PlaneSurface> surface;
  std::optional<ClusterOnTrackCalibrator> calibrator;

  void build() {
    Transform3 transform = Transform3::Identity();
    transform.linear().col(0) = Vector3::UnitY();
    transform.linear().col(1) = Vector3::UnitZ();
    transform.linear().col(2) = Vector3::UnitX();
    surface = Surface::makeShared<PlaneSurface>(
        transform,
        std::make_shared<RectangleBounds>(kBinOffset * pitch, 10_mm));
    surface->assignGeometryId(geoId);

    const auto axisX = IAxis::createEquidistant(
        AxisBoundaryType::Bound, -kBinOffset * pitch, kBinOffset * pitch,
        2 * kBinOffset, AxisDirection::AxisX);
    const auto axisY = IAxis::createEquidistant(AxisBoundaryType::Bound, -10_mm,
                                                10_mm, 1, AxisDirection::AxisY);
    DigiComponentsConfig digiCfg;
    GeometricConfig& geoCfg = digiCfg.geometricDigiConfig;
    geoCfg.indices = {eBoundLoc0};
    geoCfg.segmentation = IMultiAxis::create(*axisX, *axisY);
    geoCfg.thickness = thickness;
    geoCfg.threshold = threshold;
    geoCfg.digital = !analogue;

    ClusterOnTrackCalibrator::Config cfg;
    cfg.digitizationConfigs = DigiConfigContainer({{geoId, digiCfg}});
    cfg.chargeNoise = chargeNoise;
    cfg.minError = 0;
    calibrator.emplace(cfg);
  }

  /// The full path length through the sensor at the angle `alpha`
  double pathLength(double alpha) const { return thickness / std::cos(alpha); }

  /// Deposit the track with the position `mid` on the surface the same way
  /// as the digitization. The result maps the fired bins to their charge.
  std::map<int, double> deposit(double mid, double alpha,
                                std::mt19937& rng) const {
    const Vector4 position{0., mid, 0., 0.};
    const Vector4 momentum{std::cos(alpha), std::sin(alpha), 0., 1.};
    const ActsFatras::Hit hit(geoId, {}, position, momentum, momentum);
    const auto& geoCfg = calibrator->config()
                             .digitizationConfigs.find(geoId)
                             ->geometricDigiConfig;
    const auto channels = ActsFatras::Channelizer{}.channelize(
        hit, *surface, GeometryContext::dangerouslyDefaultConstruct(),
        Vector3::Zero(), *geoCfg.segmentation, thickness);
    BOOST_REQUIRE(channels.ok());
    std::normal_distribution<double> gauss(0., 1.);
    std::map<int, double> fired;
    for (const auto& channel : *channels) {
      double charge = channel.activation;
      if (chargeNoise > 0) {
        charge = std::max(0., charge + chargeNoise * gauss(rng));
      }
      if (charge > threshold) {
        fired[static_cast<int>(channel.bin[0]) - kBinOffset] = charge;
      }
    }
    return fired;
  }

  /// Calibrate a cluster with the predicted angle `alpha` and its error
  std::pair<double, double> calibrate(const std::map<int, double>& fired,
                                      double alpha,
                                      double alphaError = 0) const {
    Cluster cluster;
    cluster.geometryId = geoId;
    for (const auto& [bin, charge] : fired) {
      const ActsFatras::Segmentizer::Bin2D cellBin = {
          static_cast<unsigned int>(bin + kBinOffset), 0};
      cluster.channels.push_back({cellBin, {}, charge});
    }
    ClusterContainer clusters = {cluster};

    const double centre =
        0.5 * (fired.begin()->first + fired.rbegin()->first + 1) * pitch;
    MeasurementContainer measurements;
    measurements.emplaceMeasurement<1>(geoId, std::array{eBoundLoc0},
                                       Vector<1>{centre},
                                       SquareMatrix<1>{pitch * pitch / 12});

    VectorMultiTrajectory trajectory;
    auto trackState = trajectory.makeTrackState(TrackStatePropMask::Predicted |
                                                TrackStatePropMask::Calibrated);
    BoundVector parameters = BoundVector::Zero();
    parameters[eBoundPhi] = alpha;
    parameters[eBoundTheta] = 0.5 * std::numbers::pi;
    parameters[eBoundQOverP] = 1 / 1_GeV;
    BoundMatrix covariance = BoundMatrix::Identity();
    covariance(eBoundPhi, eBoundPhi) = alphaError * alphaError;
    covariance(eBoundTheta, eBoundTheta) = 0;
    trackState.predicted() = parameters;
    trackState.predictedCovariance() = covariance;
    trackState.setReferenceSurface(surface);

    const SourceLink sourceLink{IndexSourceLink(geoId, 0)};
    calibrator->calibrate(measurements, &clusters,
                          GeometryContext::dangerouslyDefaultConstruct(),
                          CalibrationContext{}, sourceLink, trackState);
    return {trackState.effectiveCalibrated()[0],
            trackState.effectiveCalibratedCovariance()(0, 0)};
  }
};

/// Fired bins from `first` to `last` with a full charge
std::map<int, double> makeBins(int first, int last) {
  std::map<int, double> fired;
  for (int bin = first; bin <= last; ++bin) {
    fired[bin] = 1;
  }
  return fired;
}

/// The RMS of the pulls of `n` tracks, uniform within one bin. Without
/// noise and without an angle error, the truth must be in the digital box.
double pullRms(const Sensor& sensor, double alpha, double alphaError,
               std::size_t n) {
  std::mt19937 rng(1234);
  std::uniform_real_distribution<double> uniform(0., sensor.pitch);
  std::normal_distribution<double> gauss(0., 1.);
  double sumPull2 = 0;
  std::size_t count = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double mid = uniform(rng);
    const auto fired = sensor.deposit(mid, alpha, rng);
    if (fired.empty()) {
      continue;
    }
    const double predicted =
        alphaError > 0 ? alpha + alphaError * gauss(rng) : alpha;
    const auto [position, variance] =
        sensor.calibrate(fired, predicted, alphaError);
    if (variance == 0) {
      // Analogue readout without noise is exact, unless the threshold
      // removed the cells after both edge bins
      CHECK_CLOSE_ABS(position, mid, 1e-9);
      continue;
    }
    if (!sensor.analogue && sensor.chargeNoise == 0 && alphaError == 0) {
      BOOST_CHECK_LE(std::abs(mid - position), std::sqrt(3 * variance) + 1e-12);
    }
    sumPull2 += std::pow(mid - position, 2) / variance;
    ++count;
  }
  return count > 0 ? std::sqrt(sumPull2 / count) : 1.;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(DigitizationSuite)

// Pitch 50 um, thickness 150 um, 10 degrees: the footprint is 0.529 pitches.
BOOST_AUTO_TEST_CASE(ClusterOnTrackDigitalBoxWidths) {
  Sensor sensor;
  sensor.build();
  const double alpha = 10_degree;
  const double phi = sensor.thickness * std::tan(alpha) / sensor.pitch;

  const auto [single, singleVar] = sensor.calibrate(makeBins(2, 2), alpha);
  CHECK_CLOSE_ABS(single, 2.5 * sensor.pitch, 1e-9);
  CHECK_CLOSE_REL(std::sqrt(singleVar),
                  (1 - phi) * sensor.pitch / std::sqrt(12), 1e-9);

  const auto [pair, pairVar] = sensor.calibrate(makeBins(2, 3), alpha);
  CHECK_CLOSE_ABS(pair, 3 * sensor.pitch, 1e-9);
  CHECK_CLOSE_REL(std::sqrt(pairVar), phi * sensor.pitch / std::sqrt(12), 1e-9);

  // Three bins cannot come from 0.529 pitches: uniform over the cluster
  const auto [triple, tripleVar] = sensor.calibrate(makeBins(2, 4), alpha);
  CHECK_CLOSE_ABS(triple, 3.5 * sensor.pitch, 1e-9);
  CHECK_CLOSE_REL(tripleVar, std::pow(3 * sensor.pitch, 2) / 12, 1e-9);
}

// The threshold shrinks the footprint by the threshold fraction at both ends
BOOST_AUTO_TEST_CASE(ClusterOnTrackDigitalThreshold) {
  Sensor sensor;
  const double alpha = 10_degree;
  sensor.threshold = 0.08 * sensor.pathLength(alpha);
  sensor.build();
  const double phi =
      sensor.thickness * std::tan(alpha) * (1 - 2 * 0.08) / sensor.pitch;

  const auto [single, singleVar] = sensor.calibrate(makeBins(0, 0), alpha);
  CHECK_CLOSE_REL(std::sqrt(singleVar),
                  (1 - phi) * sensor.pitch / std::sqrt(12), 1e-9);
  const auto [pair, pairVar] = sensor.calibrate(makeBins(0, 1), alpha);
  CHECK_CLOSE_REL(std::sqrt(pairVar), phi * sensor.pitch / std::sqrt(12), 1e-9);
}

// Without noise, the charge of the edge bins gives the exact position
BOOST_AUTO_TEST_CASE(ClusterOnTrackAnalogueExact) {
  Sensor sensor;
  sensor.analogue = true;
  sensor.build();
  const double alpha = 25_degree;
  std::mt19937 rng(42);
  const double mid = 3.3 * sensor.pitch;
  const auto [position, variance] =
      sensor.calibrate(sensor.deposit(mid, alpha, rng), alpha);
  CHECK_CLOSE_ABS(position, mid, 1e-9);
  CHECK_SMALL(variance, 1e-12);
}

BOOST_AUTO_TEST_CASE(ClusterOnTrackHighThreshold) {
  Sensor sensor;
  const double alpha = 10_degree;
  sensor.threshold = 0.75 * sensor.pathLength(alpha);
  sensor.build();
  const auto [position, variance] = sensor.calibrate(makeBins(0, 0), alpha);
  const double width = sensor.pitch - 0.5 * sensor.thickness * std::tan(alpha);
  CHECK_CLOSE_ABS(position, 0.5 * sensor.pitch, 1e-12);
  CHECK_CLOSE_REL(variance, width * width / 12, 1e-9);
}

BOOST_DATA_TEST_CASE(ClusterOnTrackSensorEdge,
                     boost::unit_test::data::make({-1., 1.}), side) {
  Sensor sensor;
  sensor.build();
  const double alpha = 10_degree;
  const double length = sensor.thickness * std::tan(alpha);
  std::mt19937 rng(42);
  const double mid = side * (kBinOffset - 0.01) * sensor.pitch;
  const auto [position, variance] =
      sensor.calibrate(sensor.deposit(mid, alpha, rng), alpha);
  const double low = (kBinOffset - 1) * sensor.pitch + 0.5 * length;
  const double high = kBinOffset * sensor.pitch;
  CHECK_CLOSE_ABS(position, side * 0.5 * (low + high), 1e-12);
  CHECK_CLOSE_REL(variance, std::pow(high - low, 2) / 12, 1e-9);
  BOOST_CHECK_LE(std::abs(position - mid), std::sqrt(3 * variance));
}

BOOST_AUTO_TEST_CASE(ClusterOnTrackInvalidConfig) {
  for (double invalid : {-1., std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
    ClusterOnTrackCalibrator::Config config;
    config.chargeNoise = invalid;
    BOOST_CHECK_THROW(ClusterOnTrackCalibrator{config}, std::invalid_argument);
    config.chargeNoise = 0;
    config.minError = invalid;
    BOOST_CHECK_THROW(ClusterOnTrackCalibrator{config}, std::invalid_argument);
  }
}

BOOST_AUTO_TEST_CASE(ClusterOnTrackInconsistentAnalogueCharge) {
  Sensor sensor;
  sensor.analogue = true;
  sensor.build();
  const double alpha = 25_degree;
  // Each edge is partial, but together their charges exceed the full path.
  const double charge = 0.6 * sensor.pathLength(alpha);
  const auto [position, variance] =
      sensor.calibrate({{0, charge}, {1, charge}}, alpha);
  CHECK_CLOSE_ABS(position, sensor.pitch, 1e-12);
  CHECK_CLOSE_REL(variance, std::pow(2 * sensor.pitch, 2) / 12, 1e-9);
}

BOOST_AUTO_TEST_CASE(ClusterOnTrackUncertainIncompatibleDirection) {
  Sensor sensor;
  sensor.chargeNoise = 1_um;
  sensor.build();
  // Normal incidence cannot give two bins, but the uncertain direction can.
  const auto [position, variance] =
      sensor.calibrate(makeBins(0, 1), 0., 5_degree);
  CHECK_CLOSE_ABS(position, sensor.pitch, 1e-9);
  BOOST_CHECK_LT(variance, sensor.pitch * sensor.pitch / 12);
}

// Without track parameters, the measurement passes through
BOOST_AUTO_TEST_CASE(ClusterOnTrackWithoutParameters) {
  Sensor sensor;
  sensor.build();
  Cluster cluster;
  cluster.channels.push_back({{kBinOffset, 0}, {}, 1});
  ClusterContainer clusters = {cluster};
  MeasurementContainer measurements;
  measurements.emplaceMeasurement<1>(sensor.geoId, std::array{eBoundLoc0},
                                     Vector<1>{0.1}, SquareMatrix<1>{0.2});
  VectorMultiTrajectory trajectory;
  auto trackState = trajectory.makeTrackState(TrackStatePropMask::Calibrated);
  trackState.setReferenceSurface(sensor.surface);
  const SourceLink sourceLink{IndexSourceLink(sensor.geoId, 0)};
  sensor.calibrator->calibrate(measurements, &clusters,
                               GeometryContext::dangerouslyDefaultConstruct(),
                               CalibrationContext{}, sourceLink, trackState);
  BOOST_CHECK_EQUAL(trackState.effectiveCalibrated()[0], 0.1);
  BOOST_CHECK_EQUAL(trackState.effectiveCalibratedCovariance()(0, 0), 0.2);
}

// The pulls have unit width for every readout
BOOST_DATA_TEST_CASE(ClusterOnTrackPulls,
                     boost::unit_test::data::make({0., 5., 10., 20., 40.}) *
                         boost::unit_test::data::make({false, true}) *
                         boost::unit_test::data::make({0., 5., 10.}),
                     angle, analogue, noise) {
  Sensor sensor;
  sensor.analogue = analogue;
  sensor.chargeNoise = noise * 1_um;
  sensor.threshold = 0.05 * sensor.pathLength(angle * 1_degree);
  sensor.build();
  CHECK_CLOSE_REL(pullRms(sensor, angle * 1_degree, 0, 20000), 1., 0.05);
}

// The error of the predicted angle is part of the calibrated error
BOOST_DATA_TEST_CASE(ClusterOnTrackAngleError,
                     boost::unit_test::data::make({10., 30.}) *
                         boost::unit_test::data::make({false, true}),
                     angle, analogue) {
  Sensor sensor;
  sensor.analogue = analogue;
  sensor.chargeNoise = 5_um;
  sensor.threshold = 0.05 * sensor.pathLength(angle * 1_degree);
  sensor.build();
  CHECK_CLOSE_REL(pullRms(sensor, angle * 1_degree, 2_degree, 20000), 1., 0.05);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
