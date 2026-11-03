// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/MbfSmoother.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cstddef>
#include <numbers>

namespace {

using namespace Acts;

using ParametersVector = Acts::BoundVector;
using CovarianceMatrix = Acts::BoundMatrix;
using Jacobian = Acts::BoundMatrix;

const Acts::GeometryContext tgContext =
    Acts::GeometryContext::dangerouslyDefaultConstruct();

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

BOOST_AUTO_TEST_CASE(Smooth) {
  VectorMultiTrajectory traj;

  std::array<BoundIndices, 2> projector{eBoundLoc0, eBoundLoc1};

  // Make dummy track parameter
  CovarianceMatrix covTrk;
  covTrk.setIdentity();
  covTrk.diagonal() << 0.08, 0.3, 1, 1, 1, 1;

  std::size_t ts_idx = traj.addTrackState(TrackStatePropMask::All);
  auto ts = traj.getTrackState(ts_idx);
  ts.typeFlags().setIsMeasurement();

  ts.predicted() << 0.3, 0.5, std::numbers::pi / 2., 1., 1 / 100., 0.;
  ts.predictedCovariance() = covTrk;

  ts.allocateCalibrated(2);
  ts.calibrated<2>() << 0.351, 0.473;
  ts.calibratedCovariance<2>() << 1e+8, 0., 0., 1e+8;
  ts.setProjectorSubspaceIndices(projector);

  ts.filtered() << 0.301, 0.503, std::numbers::pi / 2., 1., 1 / 100., 0.;
  ts.filteredCovariance() = covTrk;
  ts.pathLength() = 1.;
  ts.jacobian().setIdentity();

  ts_idx = traj.addTrackState(TrackStatePropMask::All, ts_idx);
  ts = traj.getTrackState(ts_idx);
  ts.typeFlags().setIsMeasurement();

  ts.predicted() << 0.2, 0.5, std::numbers::pi / 2., 1., 1 / 100., 0.;
  ts.predictedCovariance() = covTrk;

  ts.allocateCalibrated(2);
  ts.calibrated<2>() << 0.351, 0.473;
  ts.calibratedCovariance<2>() << 1e+8, 0., 0., 1e+8;
  ts.setProjectorSubspaceIndices(projector);

  ts.filtered() << 0.27, 0.53, std::numbers::pi / 2., 1., 1 / 100., 0.;
  ts.filteredCovariance() = covTrk;
  ts.pathLength() = 2.;
  ts.jacobian().setIdentity();

  ts_idx = traj.addTrackState(TrackStatePropMask::All, ts_idx);
  ts = traj.getTrackState(ts_idx);
  ts.typeFlags().setIsMeasurement();

  ts.predicted() << 0.35, 0.49, std::numbers::pi / 2., 1., 1 / 100., 0.;
  ts.predictedCovariance() = covTrk;

  ts.allocateCalibrated(2);
  ts.calibrated<2>() << 0.351, 0.473;
  ts.calibratedCovariance<2>() << 1e+8, 0., 0., 1e+8;
  ts.setProjectorSubspaceIndices(projector);

  ts.filtered() << 0.33, 0.43, std::numbers::pi / 2., 1., 1 / 100., 0.;
  ts.filteredCovariance() = covTrk;
  ts.pathLength() = 3.;
  ts.jacobian().setIdentity();

  // "smooth" these three track states
  BOOST_CHECK(MbfSmoother()(tgContext, traj, ts_idx).ok());

  // Regression tests, only tests very basic correctness of the math, but tests
  // for regressions in the result.

  auto ts1 = traj.getTrackState(0);
  BOOST_CHECK(ts1.hasSmoothed());
  BOOST_CHECK_NE(ts1.filtered(), ts1.smoothed());

  double tol = 1e-6;

  ParametersVector expPars;
  expPars << 0.301, 0.503, 1.5707963, 1.0, 0.01, 0.0;
  CovarianceMatrix expCov;
  expCov.setIdentity();
  expCov.diagonal() << 0.08, 0.3, 1.0, 1.0, 1.0, 1.0;
  CHECK_CLOSE_ABS(ts1.smoothed(), expPars, tol);
  CHECK_CLOSE_ABS(ts1.smoothedCovariance(), expCov, tol);

  auto ts2 = traj.getTrackState(1);
  BOOST_CHECK(ts2.hasSmoothed());
  BOOST_CHECK_NE(ts2.filtered(), ts2.smoothed());

  expPars << 0.27, 0.53, 1.5707963, 1.0, 0.01, 0.0;
  CHECK_CLOSE_ABS(ts2.smoothed(), expPars, tol);
  CHECK_CLOSE_ABS(ts2.smoothedCovariance(), expCov, tol);

  auto ts3 = traj.getTrackState(2);
  BOOST_CHECK(ts3.hasSmoothed());
  // last one, smoothed == filtered
  BOOST_CHECK_EQUAL(ts3.filtered(), ts3.smoothed());

  expPars << 0.33, 0.43, 1.5707963, 1.0, 0.01, 0.0;
  CHECK_CLOSE_ABS(ts3.smoothed(), expPars, tol);
  CHECK_CLOSE_ABS(ts3.smoothedCovariance(), expCov, tol);
}

BOOST_AUTO_TEST_CASE(MbfAgreesWithRtsAcrossCoordinateBranches) {
  for (bool phiSeam : {false, true}) {
    const auto makeTrajectory = [&] {
      VectorMultiTrajectory trajectory;
      const auto surface =
          Surface::makeShared<PlaneSurface>(Transform3::Identity());
      const auto first = trajectory.addTrackState(TrackStatePropMask::All);
      auto source = trajectory.getTrackState(first);
      source.setReferenceSurface(surface);
      source.filtered() =
          BoundVector(0., 0., phiSeam ? std::numbers::pi - 0.1 : 0.2,
                      phiSeam ? 1. : 0.125, 0.1, 0.);
      source.filteredCovariance().setIdentity();
      if (!phiSeam) {
        source.filteredCovariance()(eBoundLoc0, eBoundDir1) = 0.5;
        source.filteredCovariance()(eBoundDir1, eBoundLoc0) = 0.5;
      }
      source.predicted() = source.filtered();
      source.predictedCovariance() = source.filteredCovariance();
      auto previous = first;
      for (int step = 0; step < 3; ++step) {
        const auto index =
            trajectory.addTrackState(TrackStatePropMask::All, previous);
        auto ts = trajectory.getTrackState(index);
        const auto prev = trajectory.getTrackState(previous);
        ts.setReferenceSurface(surface);
        ts.predicted() = prev.filtered();
        ts.predictedCovariance() = prev.filteredCovariance();
        ts.jacobian().setIdentity();
        if (step == 1) {
          // Exercise backward information transport through a hole as well.
          ts.filtered() = ts.predicted();
          ts.filteredCovariance() = ts.predictedCovariance();
          previous = index;
          continue;
        }
        ts.typeFlags().setIsMeasurement();
        ts.allocateCalibrated(1);
        const BoundIndices measured =
            step == 0 ? (phiSeam ? eBoundDir0 : eBoundLoc0) : eBoundDir1;
        ts.setProjectorSubspaceIndices(std::array{measured});
        ts.calibrated<1>()[0] = step == 0
                                    ? (phiSeam ? -std::numbers::pi + 0.1 : -1.)
                                    : (phiSeam ? 1.1 : 0.25);
        ts.calibratedCovariance<1>()(0, 0) = 1.;
        BOOST_REQUIRE(GainMatrixUpdater()
                          .operator()<VectorMultiTrajectory>(tgContext, ts)
                          .ok());
        previous = index;
      }
      return trajectory;
    };
    auto mbf = makeTrajectory();
    auto rts = makeTrajectory();
    BOOST_REQUIRE(MbfSmoother()(tgContext, mbf, 3).ok());
    BOOST_REQUIRE(GainMatrixSmoother()(tgContext, rts, 3).ok());
    for (std::size_t index = 0; index < 4; ++index) {
      CHECK_CLOSE_ABS(mbf.getTrackState(index).smoothed(),
                      rts.getTrackState(index).smoothed(), 1e-12);
      CHECK_CLOSE_ABS(mbf.getTrackState(index).smoothedCovariance(),
                      rts.getTrackState(index).smoothedCovariance(), 1e-12);
    }
    if (!phiSeam) {
      const auto source = mbf.getTrackState(0);
      CHECK_CLOSE_ABS(source.smoothed()[eBoundLoc0], -31. / 60., 1e-12);
      CHECK_CLOSE_ABS(source.smoothed()[eBoundDir1], 11. / 60., 1e-12);
      CHECK_CLOSE_ABS(source.smoothedCovariance()(eBoundLoc0, eBoundDir1),
                      -2. / 15., 1e-12);
    }
  }
}

BOOST_AUTO_TEST_CASE(MbfReportsFailuresAndStopsTraversal) {
  for (bool pole : {false, true}) {
    VectorMultiTrajectory trajectory;
    const auto first = trajectory.addTrackState(TrackStatePropMask::All);
    const auto last = trajectory.addTrackState(TrackStatePropMask::All, first);
    for (auto index : {first, last}) {
      auto ts = trajectory.getTrackState(index);
      ts.filtered() = BoundVector(0., 0., 0.2, 1., 0.1, 0.);
      ts.filteredCovariance().setIdentity();
      ts.predicted() = ts.filtered();
      ts.predictedCovariance().setZero();
      ts.jacobian().setIdentity();
      ts.smoothed().setConstant(42.);
      ts.smoothedCovariance().setConstant(43.);
    }
    auto ts = trajectory.getTrackState(last);
    if (pole) {
      ts.filtered()[eBoundDir1] = 0.;
    } else {
      ts.typeFlags().setIsMeasurement();
      ts.allocateCalibrated(1);
      ts.setProjectorSubspaceIndices(std::array{eBoundLoc0});
      ts.calibrated<1>().setZero();
      ts.calibratedCovariance<1>().setZero();
    }
    const auto result = MbfSmoother()(tgContext, trajectory, last);
    BOOST_REQUIRE(!result.ok());
    if (pole) {
      BOOST_CHECK(result.error() ==
                  BoundParametrizationError::UndefinedCanonicalization);
      BOOST_CHECK(ts.smoothed().isConstant(42.));
      BOOST_CHECK(ts.smoothedCovariance().isConstant(43.));
    } else {
      BOOST_CHECK(result.error() == KalmanFitterError::SmoothFailed);
    }
    BOOST_CHECK(trajectory.getTrackState(first).smoothed().isConstant(42.));
    BOOST_CHECK(
        trajectory.getTrackState(first).smoothedCovariance().isConstant(43.));
  }
}

BOOST_AUTO_TEST_CASE(MbfDoesNotIntroduceNearAxisCutoff) {
  for (bool attachSurface : {false, true}) {
    for (double theta : {1e-14, std::numbers::pi - 1e-14}) {
      VectorMultiTrajectory trajectory;
      const auto index = trajectory.addTrackState(TrackStatePropMask::All);
      auto ts = trajectory.getTrackState(index);
      if (attachSurface) {
        ts.setReferenceSurface(
            Surface::makeShared<PlaneSurface>(Transform3::Identity()));
      }
      ts.filtered() = BoundVector(0., 0., 0.2, theta, 0.1, 0.);
      ts.filteredCovariance().setIdentity();
      BOOST_REQUIRE(MbfSmoother()(tgContext, trajectory, index).ok());
      BOOST_CHECK_GT(ts.smoothed()[eBoundDir1], 0.);
      BOOST_CHECK_LT(ts.smoothed()[eBoundDir1], std::numbers::pi);
      CHECK_CLOSE_ABS(ts.smoothed(), ts.filtered(), 1e-14);
      CHECK_CLOSE_ABS(ts.smoothedCovariance(), ts.filteredCovariance(), 1e-14);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
