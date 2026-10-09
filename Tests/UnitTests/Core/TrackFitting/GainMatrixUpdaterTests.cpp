// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <utility>

namespace {

using namespace Acts;
using namespace Acts::detail::Test;

using ParametersVector = Acts::BoundVector;
using CovarianceMatrix = Acts::BoundMatrix;
using Jacobian = Acts::BoundMatrix;

constexpr double tol = 1e-6;
const Acts::GeometryContext tgContext =
    Acts::GeometryContext::dangerouslyDefaultConstruct();

// Retain the original dense-projector calculation as an independent reference.
template <std::size_t N>
bool denseProjectorUpdate(VectorMultiTrajectory::TrackStateProxy state,
                          bool joseph) {
  const auto calibrated = state.template calibrated<N>();
  const auto calibratedCovariance = state.template calibratedCovariance<N>();
  const FixedBoundSubspaceHelper<N> subspace(
      state.projectorSubspaceIndices<N>());
  const auto H = subspace.projector();
  auto filtered = state.filtered();
  auto filteredCovariance = state.filteredCovariance();
  const auto predicted = state.predicted();
  const auto predictedCovariance = state.predictedCovariance();
  const auto K =
      (predictedCovariance * H.transpose() *
       (H * predictedCovariance * H.transpose() + calibratedCovariance)
           .inverse())
          .eval();
  if (K.hasNaN()) {
    return false;
  }
  filtered = predicted + K * (calibrated - H * predicted);
  filtered = normalizeBoundParameters(filtered);
  const auto tmp = (BoundMatrix::Identity() - K * H).eval();
  if (!joseph) {
    filteredCovariance = tmp * predictedCovariance;
  } else {
    filteredCovariance = tmp * predictedCovariance * tmp.transpose() +
                         K * calibratedCovariance * K.transpose();
  }
  const Vector<N> residual = calibrated - H * filtered;
  const SquareMatrix<N> m =
      (SquareMatrix<N>::Identity() - H * K) * calibratedCovariance;
  state.chi2() = (residual.transpose() * m.inverse() * residual).value();
  return true;
}

template <std::size_t N>
void compareCovarianceProjection() {
  std::mt19937 random(1729);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  for (std::uint8_t first = 0; first < eBoundSize; ++first) {
    for (std::uint8_t second = 0; second < eBoundSize; ++second) {
      if constexpr (N == 1) {
        if (second != 0) {
          continue;
        }
      } else if (first == second) {
        continue;
      }
      std::array<std::uint8_t, N> indices{};
      indices[0] = first;
      if constexpr (N == 2) {
        indices[1] = second;
      }
      for (bool joseph : {false, true}) {
        for (unsigned int sample = 0; sample < 32; ++sample) {
          VectorMultiTrajectory trajectory;
          auto expected = trajectory.makeTrackState(TrackStatePropMask::All);
          auto actual = trajectory.makeTrackState(TrackStatePropMask::All);
          BoundMatrix a;
          for (auto& value : a.reshaped()) {
            value = uniform(random);
          }
          BoundMatrix covariance = a * a.transpose();
          covariance.diagonal().array() += 0.01;
          if (sample % 4 == 0) {
            covariance.setZero();
            covariance.diagonal().setOnes();
          }
          if (sample == 28) {
            covariance.setZero();
          } else if (sample == 29) {
            covariance(5, 5) = std::numeric_limits<double>::infinity();
          } else if (sample == 30) {
            covariance(0, 0) = std::numeric_limits<double>::quiet_NaN();
          }
          BoundVector parameters;
          for (auto& value : parameters) {
            value = uniform(random);
          }
          Vector<N> measurement;
          for (auto& value : measurement) {
            value = uniform(random);
          }
          const SquareMatrix<N> measurementCovariance =
              sample == 28 ? SquareMatrix<N>::Zero().eval()
                           : (SquareMatrix<N>::Identity() * 0.04).eval();
          for (auto state : {expected, actual}) {
            state.predicted() = parameters;
            state.predictedCovariance() = covariance;
            state.filtered().setConstant(-111.0);
            state.filteredCovariance().setConstant(-222.0);
            state.chi2() = -333.0;
            state.allocateCalibrated(N);
            state.template calibrated<N>() = measurement;
            state.template calibratedCovariance<N>() = measurementCovariance;
            state.setProjectorSubspaceIndices(indices);
          }
          const bool expectedOk = denseProjectorUpdate<N>(expected, joseph);
          const auto result =
              GainMatrixUpdater(joseph).operator()<VectorMultiTrajectory>(
                  tgContext, actual);
          BOOST_REQUIRE_EQUAL(result.ok(), expectedOk);
          BOOST_CHECK_EQUAL(actual.chi2(), expected.chi2());
          for (unsigned int row = 0; row < eBoundSize; ++row) {
            BOOST_CHECK_EQUAL(actual.filtered()[row], expected.filtered()[row]);
            for (unsigned int col = 0; col < eBoundSize; ++col) {
              BOOST_CHECK_EQUAL(actual.filteredCovariance()(row, col),
                                expected.filteredCovariance()(row, col));
            }
          }
        }
      }
    }
  }
}

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

BOOST_AUTO_TEST_CASE(CovarianceProjectionMatchesDenseProjector) {
  compareCovarianceProjection<1>();
  compareCovarianceProjection<2>();
}

BOOST_AUTO_TEST_CASE(Update) {
  // Make dummy measurement
  Vector2 measPar(-0.1, 0.45);
  SquareMatrix2 measCov = Vector2(0.04, 0.1).asDiagonal();
  auto sourceLink = TestSourceLink(eBoundLoc0, eBoundLoc1, measPar, measCov);

  // Make dummy track parameters
  ParametersVector trkPar;
  trkPar << 0.3, 0.5, std::numbers::pi / 2., 0.3 * std::numbers::pi, 0.01, 0.;
  CovarianceMatrix trkCov = CovarianceMatrix::Zero();
  trkCov.diagonal() << 0.08, 0.3, 1, 1, 1, 0;

  // Make trajectory w/ one state
  VectorMultiTrajectory traj;
  auto idx = traj.addTrackState(TrackStatePropMask::All);
  auto ts = traj.getTrackState(idx);

  // Fill the state w/ the dummy information
  ts.predicted() = trkPar;
  ts.predictedCovariance() = trkCov;
  ts.pathLength() = 0.;
  BOOST_CHECK(!ts.hasUncalibratedSourceLink());
  testSourceLinkCalibrator<VectorMultiTrajectory>(
      tgContext, CalibrationContext{}, SourceLink{std::move(sourceLink)}, ts);
  BOOST_CHECK(ts.hasUncalibratedSourceLink());

  // Check that the state has storage available
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK(ts.hasCalibrated());

  // Gain matrix update and filtered state
  BOOST_CHECK(GainMatrixUpdater()
                  .operator()<VectorMultiTrajectory>(tgContext, ts)
                  .ok());

  // Check for regression. This does NOT test if the math is correct, just that
  // the result is the same as when the test was written.

  ParametersVector expPar;
  expPar << 0.0333333, 0.4625000, 1.5707963, 0.9424778, 0.0100000, 0.0000000;
  CHECK_CLOSE_ABS(ts.filtered(), expPar, tol);

  CovarianceMatrix expCov = CovarianceMatrix::Zero();
  expCov.diagonal() << 0.0266667, 0.0750000, 1.0000000, 1.0000000, 1.0000000,
      0.0000000;
  CHECK_CLOSE_ABS(ts.filteredCovariance(), expCov, tol);

  CHECK_CLOSE_ABS(ts.chi2(), 1.33958, 1e-4);
}

BOOST_AUTO_TEST_CASE(UpdateFailed) {
  // Make dummy measurement
  Vector2 measPar(-0.1, 0.45);
  SquareMatrix2 measCov = SquareMatrix2::Zero();
  auto sourceLink = TestSourceLink(eBoundLoc0, eBoundLoc1, measPar, measCov);

  // Make dummy track parameters
  ParametersVector trkPar;
  trkPar << 0.3, 0.5, std::numbers::pi / 2., 0.3 * std::numbers::pi, 0.01, 0.;
  CovarianceMatrix trkCov = CovarianceMatrix::Zero();

  // Make trajectory w/ one state
  VectorMultiTrajectory traj;
  auto idx = traj.addTrackState(TrackStatePropMask::All);
  auto ts = traj.getTrackState(idx);

  // Fill the state w/ the dummy information
  ts.predicted() = trkPar;
  ts.predictedCovariance() = trkCov;
  ts.pathLength() = 0.;
  BOOST_CHECK(!ts.hasUncalibratedSourceLink());
  testSourceLinkCalibrator<VectorMultiTrajectory>(
      tgContext, CalibrationContext{}, SourceLink{std::move(sourceLink)}, ts);
  BOOST_CHECK(ts.hasUncalibratedSourceLink());

  // Check that the state has storage available
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK(ts.hasCalibrated());

  // Gain matrix update and filtered state
  BOOST_CHECK(GainMatrixUpdater()
                  .operator()<VectorMultiTrajectory>(tgContext, ts)
                  .error() == KalmanFitterError::UpdateFailed);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
