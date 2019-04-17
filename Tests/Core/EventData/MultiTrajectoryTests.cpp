// This file is part of the Acts project.

// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Local Trajectory Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// clang-format on

#include <iostream>

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

using std::cout;
using std::endl;

namespace Acts {
namespace Test {

  GeometryContext gctx;

  using SourceLink = MinimalSourceLink;
  using Parameters = TrackParametersBase::ParVector_t;
  using Covariance = TrackParametersBase::CovMatrix_t;

  CurvilinearParameters
  make_params()
  {
    // generate arbitrary positive, definite matrix
    Covariance rnd = Covariance::Random();
    Covariance cov = rnd.transpose() * rnd;
    return {std::make_unique<Covariance>(cov),
            Vector3D(0, 0, 1),
            Vector3D(100, 1000, 400),
            -1};
  }

  TrackState<SourceLink, CurvilinearParameters>
  make_trackstate()
  {
    return {make_params()};
  }

  // BOOST_AUTO_TEST_CASE(multitrajectory_build)
  //{
  // MultiTrajectory<SourceLink> t;

  //// construct trajectory w/ multiple components
  // auto i0 = t.addTrackState(make_trackstate());
  //// trajectory bifurcates here into multiple hypotheses
  // auto i1a = t.addTrackState(make_trackstate(), i0);
  // auto i1b = t.addTrackState(make_trackstate(), i0);
  // auto i2a = t.addTrackState(make_trackstate(), i1a);
  // auto i2b = t.addTrackState(make_trackstate(), i1b);

  //// print each trajectory component
  // auto print = [](auto p) {
  // cout << "  point " << p.index() << endl;
  // cout << "     params " << p.predicted().transpose() << endl;
  //};
  // cout << "trajectory starting at " << i2a << endl;
  // t.visitBackwards(i2a, print);
  // cout << "trajectory starting at " << i2b << endl;
  // t.visitBackwards(i2b, print);

  //// modify elements of the trajectory
  // t.applyBackwards(i2b, [](auto p) { p.predicted().setRandom(); });
  // cout << "modified trajectory starting at " << i2b << endl;
  // t.visitBackwards(i2b, print);

  //// print full/effective parameters/covariance for an example point
  // const auto& p = t.getTrackState(i1b);
  // cout << "data for point " << p.index() << endl;
  // cout << p.predicted().transpose() << endl;
  // cout << p.predictedCovariance() << endl;
  // cout << "has uncalibrated " << p.hasUncalibrated() << endl;
  //}

  BOOST_AUTO_TEST_CASE(storage_consistency)
  {
    auto plane = Surface::makeShared<PlaneSurface>(Vector3D{0., 0., 0.},
                                                   Vector3D{1., 0., 0.});
    MultiTrajectory<SourceLink> t;
    using TrackState = TrackState<SourceLink, BoundParameters>;

    ActsMatrixD<3, 3> mCov;
    mCov << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    FittableMeasurement<SourceLink> fm
        = Measurement<SourceLink, eLOC_0, eLOC_1, ePHI>{
            plane, {}, mCov, 2, 3, 4};

    TrackState ts{SourceLink{&fm}};

    // add parameters
    using ParVec_t = BoundParameters::ParVector_t;
    using CovMat_t = BoundParameters::CovMatrix_t;

    // predicted
    ParVec_t predPar;
    predPar << 1, 2, 3, 4, 5;

    CovMat_t predCov;
    predCov << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
        19, 20, 21, 22, 23, 24, 25;

    BoundParameters pred(
        gctx, std::make_unique<CovMat_t>(predCov), predPar, plane);

    ts.parameter.predicted = pred;

    // filtered
    ParVec_t filtPar;
    filtPar << 6, 7, 8, 9, 10;

    CovMat_t filtCov;
    filtCov << 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
        42, 43, 44, 45, 46, 47, 48, 49, 50;

    BoundParameters filt(
        gctx, std::make_unique<CovMat_t>(filtCov), filtPar, plane);

    ts.parameter.filtered = filt;

    // smoothed
    ParVec_t smotPar;
    smotPar << 11, 12, 13, 14, 15;

    CovMat_t smotCov;
    smotCov << 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
        67, 68, 69, 70, 71, 72, 73, 74, 75;

    BoundParameters smot(
        gctx, std::make_unique<CovMat_t>(smotCov), smotPar, plane);

    ts.parameter.smoothed = smot;

    // "calibrate"
    ts.measurement.calibrated = **ts.measurement.uncalibrated;

    // now put it into the collection
    t.addTrackState(ts);

    // now investigate the proxy
    auto tsProxy = t.getTrackState(0);

    // std::cout << predPar.transpose() << std::endl;
    // std::cout << pred.parameters().transpose() << std::endl;
    // std::cout << (*ts.parameter.predicted).parameters().transpose()
    //<< std::endl;
    // std::cout << tsProxy.predicted().transpose() << std::endl;

    CHECK_CLOSE_ABS(pred.parameters(), tsProxy.predicted(), 1e-9);
    CHECK_CLOSE_ABS(pred.covariance(), tsProxy.predictedCovariance(), 1e-9);

    CHECK_CLOSE_ABS(filt.parameters(), tsProxy.filtered(), 1e-9);
    CHECK_CLOSE_ABS(filt.covariance(), tsProxy.filteredCovariance(), 1e-9);

    CHECK_CLOSE_ABS(smot.parameters(), tsProxy.smoothed(), 1e-9);
    CHECK_CLOSE_ABS(smot.covariance(), tsProxy.smoothedCovariance(), 1e-9);

    std::cout << tsProxy.smoothedCovariance() << std::endl;
  }

}  // namespace Test
}  // namespace Acts
