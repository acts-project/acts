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

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <iostream>

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
                                                   Vector3D{0., 0., 1.});
    MultiTrajectory<SourceLink> t;
    using TrackState = TrackState<SourceLink, BoundParameters>;

    ActsMatrixD<3, 3> mCov;
    mCov << 1, 2, 3, 4, 5, 6, 7, 8, 9;

    ActsVectorD<3> mPar;
    mPar << 2, 3, 4;
    Measurement<SourceLink, eLOC_0, eLOC_1, eQOP> meas{
        plane, {}, mCov, mPar[0], mPar[1], mPar[2]};

    FittableMeasurement<SourceLink> fm = meas;

    SourceLink sl{&fm};
    TrackState ts{sl};

    // add parameters
    using ParVec_t = BoundParameters::ParVector_t;
    using CovMat_t = BoundParameters::CovMatrix_t;

    // predicted
    ParVec_t predPar;
    predPar << 1, 2, M_PI / 4., M_PI / 2., 5;

    CovMat_t predCov;
    predCov << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
        19, 20, 21, 22, 23, 24, 25;

    BoundParameters pred(
        gctx, std::make_unique<CovMat_t>(predCov), predPar, plane);

    ts.parameter.predicted = pred;

    // filtered
    ParVec_t filtPar;
    filtPar << 6, 7, M_PI / 4., M_PI / 2., 10;

    CovMat_t filtCov;
    filtCov << 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
        42, 43, 44, 45, 46, 47, 48, 49, 50;

    BoundParameters filt(
        gctx, std::make_unique<CovMat_t>(filtCov), filtPar, plane);

    ts.parameter.filtered = filt;

    // smoothed
    ParVec_t smotPar;
    smotPar << 11, 12, M_PI / 4., M_PI / 2., 15;

    CovMat_t smotCov;
    smotCov << 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
        67, 68, 69, 70, 71, 72, 73, 74, 75;

    BoundParameters smot(
        gctx, std::make_unique<CovMat_t>(smotCov), smotPar, plane);

    ts.parameter.smoothed = smot;

    // "calibrate", keep original source link (stack address)
    ts.measurement.calibrated
        = decltype(meas){meas.referenceSurface().getSharedPtr(),
                         sl,
                         meas.covariance(),
                         meas.parameters()[0],
                         meas.parameters()[1],
                         meas.parameters()[2]};

    // make jacobian
    CovMat_t jac;
    jac << 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 89, 90, 91, 92, 93,
        94, 95, 96, 97, 98, 99, 100, 101;

    ts.parameter.jacobian = jac;

    // now put it into the collection
    t.addTrackState(ts);

    // now investigate the proxy
    auto tsProxy = t.getTrackState(0);

    // parameters
    BOOST_CHECK(tsProxy.hasPredicted());
    BOOST_CHECK_EQUAL(predPar, tsProxy.predicted());
    BOOST_CHECK_EQUAL(predCov, tsProxy.predictedCovariance());

    BOOST_CHECK(tsProxy.hasFiltered());
    BOOST_CHECK_EQUAL(filtPar, tsProxy.filtered());
    BOOST_CHECK_EQUAL(filtCov, tsProxy.filteredCovariance());

    BOOST_CHECK(tsProxy.hasSmoothed());
    BOOST_CHECK_EQUAL(smotPar, tsProxy.smoothed());
    BOOST_CHECK_EQUAL(smotCov, tsProxy.smoothedCovariance());

    BOOST_CHECK_EQUAL(&tsProxy.referenceSurface(), &ts.referenceSurface());

    BOOST_CHECK(tsProxy.hasJacobian());
    BOOST_CHECK_EQUAL(tsProxy.jacobian(), jac);

    BOOST_CHECK(tsProxy.hasProjector());
    // projector with dynamic rows
    // should be exactly equal
    BOOST_CHECK_EQUAL(tsProxy.effectiveProjector(), meas.projector());

    CovMat_t proj;
    proj << 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0;

    // full projector, should be exactly equal
    BOOST_CHECK_EQUAL(tsProxy.projector(), proj);

    // measurement properties
    BOOST_CHECK(tsProxy.hasCalibrated());
    BOOST_CHECK_EQUAL(meas.parameters(), tsProxy.effectiveCalibrated());
    ParVec_t mParFull;
    mParFull.setZero();
    mParFull.head(3) = mPar;
    BOOST_CHECK_EQUAL(mParFull, tsProxy.calibrated());

    BOOST_CHECK_EQUAL(meas.covariance(),
                      tsProxy.effectiveCalibratedCovariance());

    CovMat_t mCovFull;
    mCovFull.setZero();
    mCovFull.topLeftCorner(3, 3) = mCov;
    BOOST_CHECK_EQUAL(mCovFull, tsProxy.calibratedCovariance());

    // calibrated links to original measurement
    BOOST_CHECK_EQUAL(sl, tsProxy.calibratedSourceLink());

    // uncalibrated **is** a SourceLink
    BOOST_CHECK(tsProxy.hasUncalibrated());
    BOOST_CHECK_EQUAL(sl, tsProxy.uncalibrated());
  }

}  // namespace Test

}  // namespace Acts
