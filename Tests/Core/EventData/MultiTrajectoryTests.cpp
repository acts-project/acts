// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
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
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <iostream>

using std::cout;
using std::endl;

namespace Acts {
namespace Test {

GeometryContext gctx;

using SourceLink = MinimalSourceLink;
using Parameters = BoundVector;
using Covariance = BoundSymMatrix;

CurvilinearParameters make_params() {
  // generate arbitrary positive, definite matrix
  Covariance rnd = Covariance::Random();
  Covariance cov = rnd.transpose() * rnd;
  return {cov, Vector3D(0, 0, 1), Vector3D(100, 1000, 400), -1, 0};
}

TrackState<SourceLink, CurvilinearParameters> make_rand_trackstate() {
  return {make_params()};
}

using ParVec_t = BoundParameters::ParVector_t;
using CovMat_t = BoundParameters::CovMatrix_t;

// std::pair<TrackState<SourceLink, BoundParameters>,
// std::unique_ptr<FittableMeasurement<SourceLink>>>
auto make_trackstate() {
  auto plane = Surface::makeShared<PlaneSurface>(Vector3D{0., 0., 0.},
                                                 Vector3D{0., 0., 1.});
  using TrackState = TrackState<SourceLink, BoundParameters>;

  ActsMatrixD<3, 3> mCov;
  mCov << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  ActsVectorD<3> mPar;
  mPar << 2, 3, 4;
  Measurement<SourceLink, eLOC_0, eLOC_1, eQOP> meas{plane,   {},      mCov,
                                                     mPar[0], mPar[1], mPar[2]};

  auto fm = std::make_unique<FittableMeasurement<SourceLink>>(meas);

  SourceLink sl{fm.get()};
  TrackState ts{sl};

  // add parameters

  // predicted
  ParVec_t predPar;
  predPar << 1, 2, M_PI / 4., M_PI / 2., 5, 0.;

  CovMat_t predCov;
  predCov.setRandom();

  BoundParameters pred(gctx, predCov, predPar, plane);

  ts.parameter.predicted = pred;

  // filtered
  ParVec_t filtPar;
  filtPar << 6, 7, M_PI / 4., M_PI / 2., 10, 0.;

  CovMat_t filtCov;
  filtCov.setRandom();

  BoundParameters filt(gctx, filtCov, filtPar, plane);

  ts.parameter.filtered = filt;

  // smoothed
  ParVec_t smotPar;
  smotPar << 11, 12, M_PI / 4., M_PI / 2., 15, 0.;

  CovMat_t smotCov;
  smotCov.setRandom();

  BoundParameters smot(gctx, smotCov, smotPar, plane);

  ts.parameter.smoothed = smot;

  // "calibrate", keep original source link (stack address)
  ts.measurement.calibrated =
      decltype(meas){meas.referenceSurface().getSharedPtr(),
                     sl,
                     meas.covariance(),
                     meas.parameters()[0],
                     meas.parameters()[1],
                     meas.parameters()[2]};

  // make jacobian
  CovMat_t jac;
  jac.setRandom();

  ts.parameter.jacobian = jac;

  ts.parameter.chi2 = 78;
  ts.parameter.pathLength = 42;

  return std::make_tuple(ts, std::move(fm), meas);
}

BOOST_AUTO_TEST_CASE(multitrajectory_build) {
  MultiTrajectory<SourceLink> t;

  // construct trajectory w/ multiple components
  auto i0 = t.addTrackState(make_rand_trackstate());
  // trajectory bifurcates here into multiple hypotheses
  auto i1a = t.addTrackState(make_rand_trackstate(), i0);
  auto i1b = t.addTrackState(make_rand_trackstate(), i0);
  auto i2a = t.addTrackState(make_rand_trackstate(), i1a);
  auto i2b = t.addTrackState(make_rand_trackstate(), i1b);

  // print each trajectory component
  std::vector<size_t> act;
  auto collect = [&](auto p) {
    act.push_back(p.index());
    // assert absence of things
    BOOST_CHECK(!p.hasUncalibrated());
    BOOST_CHECK(!p.hasCalibrated());
    BOOST_CHECK(!p.hasFiltered());
    BOOST_CHECK(!p.hasSmoothed());
    BOOST_CHECK(!p.hasJacobian());
    BOOST_CHECK(!p.hasProjector());
  };

  std::vector<size_t> exp = {i2a, i1a, i0};
  t.visitBackwards(i2a, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());

  act.clear();
  exp = {i2b, i1b, i0};
  t.visitBackwards(i2b, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());

  act.clear();
  t.applyBackwards(i2b, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());
}

BOOST_AUTO_TEST_CASE(visit_apply_abort) {
  MultiTrajectory<SourceLink> t;

  // construct trajectory with three components
  auto i0 = t.addTrackState(make_rand_trackstate());
  auto i1 = t.addTrackState(make_rand_trackstate(), i0);
  auto i2 = t.addTrackState(make_rand_trackstate(), i1);

  size_t n = 0;
  t.applyBackwards(i2, [&](const auto&) {
    n++;
    return false;
  });
  BOOST_CHECK_EQUAL(n, 1u);

  n = 0;
  t.applyBackwards(i2, [&](const auto& ts) {
    n++;
    if (ts.index() == i1) {
      return false;
    }
    return true;
  });
  BOOST_CHECK_EQUAL(n, 2u);

  n = 0;
  t.applyBackwards(i2, [&](const auto&) {
    n++;
    return true;
  });
  BOOST_CHECK_EQUAL(n, 3u);
}

BOOST_AUTO_TEST_CASE(trackstate_add_bitmask) {
  namespace PM = TrackStatePropMask;
  auto bs1 = PM::Uncalibrated;

  BOOST_CHECK(ACTS_CHECK_BIT(bs1, PM::Uncalibrated));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs1, PM::Calibrated));

  auto bs2 = PM::Calibrated;

  BOOST_CHECK(!ACTS_CHECK_BIT(bs2, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(bs2, PM::Calibrated));

  auto bs3 = PM::Calibrated | PM::Uncalibrated;

  BOOST_CHECK(ACTS_CHECK_BIT(bs3, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(bs3, PM::Calibrated));

  BOOST_CHECK(ACTS_CHECK_BIT(PM::All, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(PM::All, PM::Calibrated));

  auto bs4 = PM::Predicted | PM::Jacobian | PM::Uncalibrated;
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Predicted));
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Jacobian));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Calibrated));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Filtered));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Smoothed));

  MultiTrajectory<SourceLink> t;

  auto ts = t.getTrackState(t.addTrackState(PM::All));
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK(ts.hasUncalibrated());
  BOOST_CHECK(ts.hasCalibrated());
  BOOST_CHECK(ts.hasProjector());
  BOOST_CHECK(ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::None));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Predicted));
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Filtered));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Smoothed));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Uncalibrated));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Calibrated));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(ts.hasCalibrated());
  BOOST_CHECK(ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Jacobian));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(ts.hasJacobian());
}

BOOST_AUTO_TEST_CASE(trackstate_proxy_cross_talk) {
  // assert expected "cross-talk" between trackstate proxies
  auto [ots, fm, meas] = make_trackstate();

  MultiTrajectory<SourceLink> t;

  t.addTrackState(ots);

  const auto& ct = t;
  auto cts = ct.getTrackState(0);
  auto ts = t.getTrackState(0);

  // assert expected value of chi2 and path length
  BOOST_CHECK_EQUAL(cts.chi2(), 78u);
  BOOST_CHECK_EQUAL(ts.chi2(), 78u);
  BOOST_CHECK_EQUAL(cts.pathLength(), 42u);
  BOOST_CHECK_EQUAL(ts.pathLength(), 42u);

  ParVec_t v;
  CovMat_t cov;

  v.setRandom();
  ts.predicted() = v;
  BOOST_CHECK_EQUAL(cts.predicted(), v);
  cov.setRandom();
  ts.predictedCovariance() = cov;
  BOOST_CHECK_EQUAL(cts.predictedCovariance(), cov);

  v.setRandom();
  ts.filtered() = v;
  BOOST_CHECK_EQUAL(cts.filtered(), v);
  cov.setRandom();
  ts.filteredCovariance() = cov;
  BOOST_CHECK_EQUAL(cts.filteredCovariance(), cov);

  v.setRandom();
  ts.smoothed() = v;
  BOOST_CHECK_EQUAL(cts.smoothed(), v);
  cov.setRandom();
  ts.smoothedCovariance() = cov;
  BOOST_CHECK_EQUAL(cts.smoothedCovariance(), cov);

  // make copy of fm
  auto fm2 = std::make_unique<FittableMeasurement<SourceLink>>(*fm);
  SourceLink sl2{fm2.get()};
  ts.uncalibrated() = sl2;
  BOOST_CHECK_EQUAL(cts.uncalibrated(), sl2);
  BOOST_CHECK_NE(cts.uncalibrated(), SourceLink{fm.get()});

  CovMat_t newMeasCov;
  newMeasCov.setRandom();
  ts.calibratedCovariance() = newMeasCov;
  BOOST_CHECK_EQUAL(cts.calibratedCovariance(), newMeasCov);

  ParVec_t newMeasPar;
  newMeasPar.setRandom();
  ts.calibrated() = newMeasPar;
  BOOST_CHECK_EQUAL(cts.calibrated(), newMeasPar);

  size_t measdim = ts.effectiveCalibrated().rows();

  ActsMatrixXd eff{measdim, measdim};
  eff.setRandom();
  ts.effectiveCalibratedCovariance() = eff;
  BOOST_CHECK_EQUAL(cts.effectiveCalibratedCovariance(), eff);
  newMeasCov.topLeftCorner(eff.rows(), eff.rows()) = eff;
  BOOST_CHECK_EQUAL(cts.calibratedCovariance(), newMeasCov);

  CovMat_t jac;
  jac.setRandom();
  ts.jacobian() = jac;
  BOOST_CHECK_EQUAL(cts.jacobian(), jac);

  ts.chi2() = 98;
  BOOST_CHECK_EQUAL(cts.chi2(), 98u);

  ts.pathLength() = 66;
  BOOST_CHECK_EQUAL(cts.pathLength(), 66u);
}

BOOST_AUTO_TEST_CASE(trackstate_reassignment) {
  auto [ots, fm, meas] = make_trackstate();

  constexpr size_t maxmeasdim = MultiTrajectory<SourceLink>::MeasurementSizeMax;

  MultiTrajectory<SourceLink> t;
  t.addTrackState(ots);

  auto ts = t.getTrackState(0);

  // assert measdim and contents of original measurement (just to be safe)
  BOOST_CHECK_EQUAL(ts.calibratedSize(), meas.size());
  BOOST_CHECK_EQUAL(ts.effectiveCalibrated(), meas.parameters());
  BOOST_CHECK_EQUAL(ts.effectiveCalibratedCovariance(), meas.covariance());
  BOOST_CHECK_EQUAL(ts.effectiveProjector(), meas.projector());

  // create new measurement
  ActsSymMatrixD<2> mCov;
  mCov.setRandom();
  ActsVectorD<2> mPar;
  mPar.setRandom();
  Measurement<SourceLink, eLOC_0, eLOC_1> m2{
      meas.referenceSurface().getSharedPtr(), {}, mCov, mPar[0], mPar[1]};

  ts.setCalibrated(m2);

  BOOST_CHECK_EQUAL(ts.calibratedSize(), 2u);
  BOOST_CHECK_EQUAL(ts.effectiveCalibrated(), mPar);
  BOOST_CHECK_EQUAL(ts.effectiveCalibratedCovariance(), mCov);
  BOOST_CHECK_EQUAL(ts.effectiveProjector(), m2.projector());

  // check if overallocated part is zeroed correctly
  ActsVectorD<maxmeasdim> mParFull;
  mParFull.setZero();
  mParFull.head(2) = mPar;
  BOOST_CHECK_EQUAL(ts.calibrated(), mParFull);

  ActsSymMatrixD<maxmeasdim> mCovFull;
  mCovFull.setZero();
  mCovFull.topLeftCorner(2, 2) = mCov;
  BOOST_CHECK_EQUAL(ts.calibratedCovariance(), mCovFull);

  ActsSymMatrixD<maxmeasdim> projFull;
  projFull.setZero();
  projFull.topLeftCorner(m2.size(), maxmeasdim) = m2.projector();
  BOOST_CHECK_EQUAL(ts.projector(), projFull);
}

BOOST_AUTO_TEST_CASE(storage_consistency) {
  MultiTrajectory<SourceLink> t;

  auto [ts, fm, om] = make_trackstate();

  // now put it into the collection
  t.addTrackState(ts);

  // now investigate the proxy
  auto tsProxy = t.getTrackState(0);

  // parameters
  BOOST_CHECK(tsProxy.hasPredicted());
  BOOST_CHECK_EQUAL(ts.parameter.predicted->parameters(), tsProxy.predicted());
  BOOST_CHECK_EQUAL(*ts.parameter.predicted->covariance(),
                    tsProxy.predictedCovariance());

  BOOST_CHECK(tsProxy.hasFiltered());
  BOOST_CHECK_EQUAL(ts.parameter.filtered->parameters(), tsProxy.filtered());
  BOOST_CHECK_EQUAL(*ts.parameter.filtered->covariance(),
                    tsProxy.filteredCovariance());

  BOOST_CHECK(tsProxy.hasSmoothed());
  BOOST_CHECK_EQUAL(ts.parameter.smoothed->parameters(), tsProxy.smoothed());
  BOOST_CHECK_EQUAL(*ts.parameter.smoothed->covariance(),
                    tsProxy.smoothedCovariance());

  BOOST_CHECK_EQUAL(&tsProxy.referenceSurface(), &ts.referenceSurface());

  BOOST_CHECK(tsProxy.hasJacobian());
  BOOST_CHECK_EQUAL(tsProxy.jacobian(), *ts.parameter.jacobian);

  BOOST_CHECK(tsProxy.hasProjector());
  std::visit(
      [&](const auto& meas) {
        BOOST_CHECK_EQUAL(tsProxy.effectiveProjector(), meas.projector());
        // measurement properties
        BOOST_CHECK(tsProxy.hasCalibrated());
        BOOST_CHECK_EQUAL(meas.parameters(), tsProxy.effectiveCalibrated());
        ParVec_t mParFull;
        mParFull.setZero();
        mParFull.head(meas.size()) = meas.parameters();
        BOOST_CHECK_EQUAL(mParFull, tsProxy.calibrated());

        BOOST_CHECK_EQUAL(meas.covariance(),
                          tsProxy.effectiveCalibratedCovariance());
        CovMat_t mCovFull;
        mCovFull.setZero();
        mCovFull.topLeftCorner(meas.size(), meas.size()) = meas.covariance();
        BOOST_CHECK_EQUAL(mCovFull, tsProxy.calibratedCovariance());

        // calibrated links to original measurement
        BOOST_CHECK_EQUAL(meas.sourceLink(), tsProxy.calibratedSourceLink());

        // uncalibrated **is** a SourceLink
        BOOST_CHECK(tsProxy.hasUncalibrated());
        BOOST_CHECK_EQUAL(meas.sourceLink(), tsProxy.uncalibrated());

        // full projector, should be exactly equal
        CovMat_t fullProj;
        fullProj.setZero();
        fullProj.topLeftCorner(
            meas.size(), MultiTrajectory<SourceLink>::MeasurementSizeMax) =
            meas.projector();
        BOOST_CHECK_EQUAL(tsProxy.projector(), fullProj);

        // projector with dynamic rows
        // should be exactly equal
        BOOST_CHECK_EQUAL(tsProxy.effectiveProjector(), meas.projector());
      },
      *ts.measurement.calibrated);
}

}  // namespace Test

}  // namespace Acts
