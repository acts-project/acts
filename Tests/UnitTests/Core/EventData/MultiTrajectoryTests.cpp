// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/context.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/GenerateParameters.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <numeric>
#include <optional>
#include <random>
#include <tuple>

namespace {

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::Test;
namespace bd = boost::unit_test::data;

using ParametersVector = BoundTrackParameters::ParametersVector;
using CovarianceMatrix = BoundTrackParameters::CovarianceMatrix;
using Jacobian = BoundMatrix;

struct TestTrackState {
  std::shared_ptr<Surface> surface;
  TestSourceLink sourceLink;
  BoundTrackParameters predicted;
  BoundTrackParameters filtered;
  BoundTrackParameters smoothed;
  Jacobian jacobian;
  double chi2;
  double pathLength;

  // Generate a random TestTrackState.
  //
  // @param rng Random number generator
  // @param size_t nMeasurement either 1 or 2
  template <typename rng_t>
  TestTrackState(rng_t& rng, size_t nMeasurements)
      : surface(Surface::makeShared<PlaneSurface>(Vector3::Zero(),
                                                  Vector3::UnitZ())),
        // set bogus parameters first since they are not default-constructible
        predicted(surface, BoundVector::Zero()),
        filtered(surface, BoundVector::Zero()),
        smoothed(surface, BoundVector::Zero()),
        jacobian(Jacobian::Identity()),
        chi2(std::chi_squared_distribution<double>(nMeasurements)(rng)),
        pathLength(
            std::uniform_real_distribution<ActsScalar>(1_mm, 10_mm)(rng)) {
    // set a random geometry identifier to uniquely identify each surface
    auto geoId =
        std::uniform_int_distribution<GeometryIdentifier::Value>()(rng);
    surface->assignGeometryId(geoId);

    // create source link w/ inline 1d or 2d measurement data
    if (nMeasurements == 1u) {
      auto [par, cov] = generateParametersCovariance<ActsScalar, 1u>(rng);
      sourceLink = TestSourceLink(eBoundLoc0, par[0], cov(0, 0), geoId);
    } else if (nMeasurements == 2u) {
      auto [par, cov] = generateParametersCovariance<ActsScalar, 2u>(rng);
      sourceLink = TestSourceLink(eBoundLoc1, eBoundQOverP, par, cov, geoId);
    } else {
      throw std::runtime_error("invalid number of measurement dimensions");
    }

    // create track parameters
    auto [trkPar, trkCov] = generateBoundParametersCovariance(rng);
    // trkPar[eBoundPhi] = 45_degree;
    // trkPar[eBoundTheta] = 90_degree;
    // trkPar[eBoundQOverP] = 5.;
    // predicted
    predicted = BoundTrackParameters(surface, trkPar, trkCov);
    // filtered, modified q/p, reduced covariance
    // trkPar[eBoundQOverP] = 10.;
    filtered = BoundTrackParameters(surface, trkPar, 0.75 * trkCov);
    // smoothed, modified q/p, further reduced covariance
    // trkPar[eBoundQOverP] = 15.;
    smoothed = BoundTrackParameters(surface, trkPar, 0.5 * trkCov);

    // propagation jacobian is identity + corrections
    for (Eigen::Index c = 0; c < jacobian.cols(); ++c) {
      for (Eigen::Index r = 0; r < jacobian.rows(); ++r) {
        jacobian(c, r) +=
            std::uniform_real_distribution<ActsScalar>(-0.125, 0.125)(rng);
      }
    }
  }
};

// Fill a TrackStateProxy with values from a TestTrackState.
//
// @param[in] pc TestTrackState with the input values
// @param[in] mask Specifies which components are used/filled
// @param[out] ts TrackStateProxy which is filled
// @param [in] nMeasurements Dimension of the measurement
template <typename track_state_t>
void fillTrackState(const TestTrackState& pc, TrackStatePropMask mask,
                    track_state_t& ts) {
  // always set the reference surface
  ts.setReferenceSurface(pc.predicted.referenceSurface().getSharedPtr());

  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Predicted)) {
    ts.predicted() = pc.predicted.parameters();
    BOOST_CHECK(pc.predicted.covariance().has_value());
    ts.predictedCovariance() = *(pc.predicted.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Filtered)) {
    ts.filtered() = pc.filtered.parameters();
    BOOST_CHECK(pc.filtered.covariance().has_value());
    ts.filteredCovariance() = *(pc.filtered.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Smoothed)) {
    ts.smoothed() = pc.smoothed.parameters();
    BOOST_CHECK(pc.smoothed.covariance().has_value());
    ts.smoothedCovariance() = *(pc.smoothed.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Jacobian)) {
    ts.jacobian() = pc.jacobian;
  }
  ts.chi2() = pc.chi2;
  ts.pathLength() = pc.pathLength;
  // source link defines the uncalibrated measurement
  ts.setUncalibrated(pc.sourceLink);
  // create calibrated measurements from source link
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Calibrated)) {
    testSourceLinkCalibrator<VectorMultiTrajectory>(Acts::GeometryContext{},
                                                    ts);
  }
}

const GeometryContext gctx;
// fixed seed for reproducible tests
std::default_random_engine rng(31415);

}  // namespace

BOOST_AUTO_TEST_SUITE(EventDataMultiTrajectory)

BOOST_AUTO_TEST_CASE(Build) {
  constexpr TrackStatePropMask kMask = TrackStatePropMask::Predicted;

  // construct trajectory w/ multiple components
  VectorMultiTrajectory t;

  auto i0 = t.addTrackState(kMask);
  // trajectory bifurcates here into multiple hypotheses
  auto i1a = t.addTrackState(kMask, i0);
  auto i1b = t.addTrackState(kMask, i0);
  auto i2a = t.addTrackState(kMask, i1a);
  auto i2b = t.addTrackState(kMask, i1b);

  // print each trajectory component
  std::vector<size_t> act;
  auto collect = [&](auto p) {
    act.push_back(p.index());
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

BOOST_AUTO_TEST_CASE(AsReadOnly) {
  // make mutable
  VectorMultiTrajectory t;
  auto i0 = t.addTrackState();

  BOOST_CHECK(!t.ReadOnly);

  {
    VectorMultiTrajectory::TrackStateProxy tsp = t.getTrackState(i0);
    static_cast<void>(tsp);
    VectorMultiTrajectory::ConstTrackStateProxy ctsp = t.getTrackState(i0);
    static_cast<void>(ctsp);
  }

  ConstVectorMultiTrajectory ct = t;

  ConstVectorMultiTrajectory ctm{std::move(t)};

  {
    static_assert(
        std::is_same_v<ConstVectorMultiTrajectory::ConstTrackStateProxy,
                       decltype(ct.getTrackState(i0))>,
        "Got mutable track state proxy");
    ConstVectorMultiTrajectory::ConstTrackStateProxy ctsp =
        ct.getTrackState(i0);
    static_cast<void>(ctsp);

    // doesn't compile:
    // ctsp.predictedCovariance().setIdentity();
  }

  // doesn't compile:
  // ct.clear();
  // ct.addTrackState();
}

BOOST_AUTO_TEST_CASE(Clear) {
  constexpr TrackStatePropMask kMask = TrackStatePropMask::Predicted;
  VectorMultiTrajectory t;
  BOOST_CHECK_EQUAL(t.size(), 0);

  auto i0 = t.addTrackState(kMask);
  // trajectory bifurcates here into multiple hypotheses
  auto i1a = t.addTrackState(kMask, i0);
  auto i1b = t.addTrackState(kMask, i0);
  t.addTrackState(kMask, i1a);
  t.addTrackState(kMask, i1b);

  BOOST_CHECK_EQUAL(t.size(), 5);
  t.clear();
  BOOST_CHECK_EQUAL(t.size(), 0);
}

BOOST_AUTO_TEST_CASE(ApplyWithAbort) {
  constexpr TrackStatePropMask kMask = TrackStatePropMask::Predicted;

  // construct trajectory with three components
  VectorMultiTrajectory t;
  auto i0 = t.addTrackState(kMask);
  auto i1 = t.addTrackState(kMask, i0);
  auto i2 = t.addTrackState(kMask, i1);

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

BOOST_AUTO_TEST_CASE(BitmaskOperators) {
  using PM = TrackStatePropMask;

  auto bs1 = PM::Predicted;

  BOOST_CHECK(ACTS_CHECK_BIT(bs1, PM::Predicted));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs1, PM::Calibrated));

  auto bs2 = PM::Calibrated;

  BOOST_CHECK(!ACTS_CHECK_BIT(bs2, PM::Predicted));
  BOOST_CHECK(ACTS_CHECK_BIT(bs2, PM::Calibrated));

  auto bs3 = PM::Calibrated | PM::Predicted;

  BOOST_CHECK(ACTS_CHECK_BIT(bs3, PM::Predicted));
  BOOST_CHECK(ACTS_CHECK_BIT(bs3, PM::Calibrated));

  BOOST_CHECK(ACTS_CHECK_BIT(PM::All, PM::Predicted));
  BOOST_CHECK(ACTS_CHECK_BIT(PM::All, PM::Calibrated));

  auto bs4 = PM::Predicted | PM::Jacobian | PM::Smoothed;
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Predicted));
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Jacobian));
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Smoothed));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Calibrated));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Filtered));

  auto cnv = [](auto a) -> std::bitset<8> {
    return static_cast<std::underlying_type<PM>::type>(a);
  };

  BOOST_CHECK(cnv(PM::All).all());    // all ones
  BOOST_CHECK(cnv(PM::None).none());  // all zeros

  // test orthogonality
  std::array<PM, 5> values{PM::Predicted, PM::Filtered, PM::Smoothed,
                           PM::Jacobian, PM::Calibrated};
  for (size_t i = 0; i < values.size(); i++) {
    for (size_t j = 0; j < values.size(); j++) {
      PM a = values[i];
      PM b = values[j];

      if (i == j) {
        BOOST_CHECK(cnv(a & b).count() == 1);
      } else {
        BOOST_CHECK(cnv(a & b).none());
      }
    }
  }

  BOOST_CHECK(cnv(PM::Predicted ^ PM::Filtered).count() == 2);
  BOOST_CHECK(cnv(PM::Predicted ^ PM::Predicted).none());
  BOOST_CHECK(~(PM::Predicted | PM::Calibrated) ==
              (PM::All ^ PM::Predicted ^ PM::Calibrated));

  PM base = PM::None;
  BOOST_CHECK(cnv(base) == 0);

  base &= PM::Filtered;
  BOOST_CHECK(cnv(base) == 0);

  base |= PM::Filtered;
  BOOST_CHECK(base == PM::Filtered);

  base |= PM::Calibrated;
  BOOST_CHECK(base == (PM::Filtered | PM::Calibrated));

  base ^= PM::All;
  BOOST_CHECK(base == ~(PM::Filtered | PM::Calibrated));
}

BOOST_AUTO_TEST_CASE(AddTrackStateWithBitMask) {
  using PM = TrackStatePropMask;
  using namespace Acts::HashedStringLiteral;

  VectorMultiTrajectory t;

  auto alwaysPresent = [](auto& ts) {
    BOOST_CHECK(ts.template has<"calibratedSourceLink"_hash>());
    BOOST_CHECK(ts.template has<"referenceSurface"_hash>());
    BOOST_CHECK(ts.template has<"measdim"_hash>());
    BOOST_CHECK(ts.template has<"chi2"_hash>());
    BOOST_CHECK(ts.template has<"pathLength"_hash>());
    BOOST_CHECK(ts.template has<"typeFlags"_hash>());
  };

  auto ts = t.getTrackState(t.addTrackState(PM::All));
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(ts.hasProjector());
  BOOST_CHECK(ts.hasJacobian());
  alwaysPresent(ts);
  ts.allocateCalibrated(5);
  BOOST_CHECK(ts.hasCalibrated());

  ts = t.getTrackState(t.addTrackState(PM::None));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());
  alwaysPresent(ts);

  ts = t.getTrackState(t.addTrackState(PM::Predicted));
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());
  alwaysPresent(ts);

  ts = t.getTrackState(t.addTrackState(PM::Filtered));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());
  alwaysPresent(ts);

  ts = t.getTrackState(t.addTrackState(PM::Smoothed));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());
  alwaysPresent(ts);

  ts = t.getTrackState(t.addTrackState(PM::Calibrated));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());
  ts.allocateCalibrated(5);
  BOOST_CHECK(ts.hasCalibrated());

  ts = t.getTrackState(t.addTrackState(PM::Jacobian));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(ts.hasJacobian());
  alwaysPresent(ts);
}

// assert expected "cross-talk" between trackstate proxies
BOOST_AUTO_TEST_CASE(TrackStateProxyCrossTalk) {
  TestTrackState pc(rng, 2u);

  // multi trajectory w/ a single, fully set, track state
  VectorMultiTrajectory traj;
  size_t index = traj.addTrackState();
  {
    auto ts = traj.getTrackState(index);
    fillTrackState(pc, TrackStatePropMask::All, ts);
  }
  // get two TrackStateProxies that reference the same data
  auto tsa = traj.getTrackState(index);
  auto tsb = traj.getTrackState(index);
  // then modify one and check that the other was modified as well
  {
    auto [par, cov] = generateBoundParametersCovariance(rng);
    tsb.predicted() = par;
    tsb.predictedCovariance() = cov;
    BOOST_CHECK_EQUAL(tsa.predicted(), par);
    BOOST_CHECK_EQUAL(tsa.predictedCovariance(), cov);
    BOOST_CHECK_EQUAL(tsb.predicted(), par);
    BOOST_CHECK_EQUAL(tsb.predictedCovariance(), cov);
  }
  {
    auto [par, cov] = generateBoundParametersCovariance(rng);
    tsb.filtered() = par;
    tsb.filteredCovariance() = cov;
    BOOST_CHECK_EQUAL(tsa.filtered(), par);
    BOOST_CHECK_EQUAL(tsa.filteredCovariance(), cov);
    BOOST_CHECK_EQUAL(tsb.filtered(), par);
    BOOST_CHECK_EQUAL(tsb.filteredCovariance(), cov);
  }
  {
    auto [par, cov] = generateBoundParametersCovariance(rng);
    tsb.smoothed() = par;
    tsb.smoothedCovariance() = cov;
    BOOST_CHECK_EQUAL(tsa.smoothed(), par);
    BOOST_CHECK_EQUAL(tsa.smoothedCovariance(), cov);
    BOOST_CHECK_EQUAL(tsb.smoothed(), par);
    BOOST_CHECK_EQUAL(tsb.smoothedCovariance(), cov);
  }
  {
    // create a new (invalid) source link
    TestSourceLink invalid;
    BOOST_CHECK_NE(&tsa.uncalibrated(), &invalid);
    BOOST_CHECK_NE(&tsb.uncalibrated(), &invalid);
    tsb.setUncalibrated(invalid);
    BOOST_CHECK_EQUAL(&tsa.uncalibrated(), &invalid);
    BOOST_CHECK_EQUAL(&tsb.uncalibrated(), &invalid);
  }
  {
    // reset measurements w/ full parameters
    auto [measPar, measCov] = generateBoundParametersCovariance(rng);
    tsb.allocateCalibrated(eBoundSize);
    tsb.calibrated<eBoundSize>() = measPar;
    tsb.calibratedCovariance<eBoundSize>() = measCov;
    BOOST_CHECK_EQUAL(tsa.calibrated<eBoundSize>(), measPar);
    BOOST_CHECK_EQUAL(tsa.calibratedCovariance<eBoundSize>(), measCov);
    BOOST_CHECK_EQUAL(tsb.calibrated<eBoundSize>(), measPar);
    BOOST_CHECK_EQUAL(tsb.calibratedCovariance<eBoundSize>(), measCov);
  }
  {
    // reset only the effective measurements
    auto [measPar, measCov] = generateBoundParametersCovariance(rng);
    size_t nMeasurements = tsb.effectiveCalibrated().rows();
    auto effPar = measPar.head(nMeasurements);
    auto effCov = measCov.topLeftCorner(nMeasurements, nMeasurements);
    tsb.allocateCalibrated(eBoundSize);
    tsb.effectiveCalibrated() = effPar;
    tsb.effectiveCalibratedCovariance() = effCov;
    BOOST_CHECK_EQUAL(tsa.effectiveCalibrated(), effPar);
    BOOST_CHECK_EQUAL(tsa.effectiveCalibratedCovariance(), effCov);
    BOOST_CHECK_EQUAL(tsb.effectiveCalibrated(), effPar);
    BOOST_CHECK_EQUAL(tsb.effectiveCalibratedCovariance(), effCov);
  }
  {
    Jacobian jac = Jacobian::Identity();
    BOOST_CHECK_NE(tsa.jacobian(), jac);
    BOOST_CHECK_NE(tsb.jacobian(), jac);
    tsb.jacobian() = jac;
    BOOST_CHECK_EQUAL(tsa.jacobian(), jac);
    BOOST_CHECK_EQUAL(tsb.jacobian(), jac);
  }
  {
    tsb.chi2() = 98.0;
    BOOST_CHECK_EQUAL(tsa.chi2(), 98.0);
    BOOST_CHECK_EQUAL(tsb.chi2(), 98.0);
  }
  {
    tsb.pathLength() = 66.0;
    BOOST_CHECK_EQUAL(tsa.pathLength(), 66.0);
    BOOST_CHECK_EQUAL(tsb.pathLength(), 66.0);
  }
}

BOOST_AUTO_TEST_CASE(TrackStateReassignment) {
  TestTrackState pc(rng, 1u);

  VectorMultiTrajectory t;
  size_t index = t.addTrackState();
  auto ts = t.getTrackState(index);
  fillTrackState(pc, TrackStatePropMask::All, ts);

  // assert contents of original measurement (just to be safe)
  BOOST_CHECK_EQUAL(ts.calibratedSize(), 1u);
  BOOST_CHECK_EQUAL(ts.effectiveCalibrated(),
                    (pc.sourceLink.parameters.head<1>()));
  BOOST_CHECK_EQUAL(ts.effectiveCalibratedCovariance(),
                    (pc.sourceLink.covariance.topLeftCorner<1, 1>()));

  // use temporary measurement to reset calibrated data
  TestTrackState ttsb(rng, 2u);
  ts.setUncalibrated(ttsb.sourceLink);
  auto meas = testSourceLinkCalibratorReturn<VectorMultiTrajectory>(gctx, ts);
  auto m2 = std::get<Measurement<BoundIndices, 2u>>(meas);

  BOOST_CHECK_EQUAL(ts.calibratedSize(), 2);
  BOOST_CHECK_EQUAL(ts.effectiveCalibrated(), m2.parameters());
  BOOST_CHECK_EQUAL(ts.effectiveCalibratedCovariance(), m2.covariance());
  BOOST_CHECK_EQUAL(ts.effectiveProjector(), m2.projector());
}

BOOST_DATA_TEST_CASE(TrackStateProxyStorage, bd::make({1u, 2u}),
                     nMeasurements) {
  TestTrackState pc(rng, nMeasurements);

  // create trajectory with a single fully-filled random track state
  VectorMultiTrajectory t;
  size_t index = t.addTrackState();
  auto ts = t.getTrackState(index);
  fillTrackState(pc, TrackStatePropMask::All, ts);

  // check that the surface is correctly set
  BOOST_CHECK_EQUAL(&ts.referenceSurface(), pc.surface.get());
  BOOST_CHECK_EQUAL(ts.referenceSurface().geometryId(),
                    pc.sourceLink.geometryId());

  // check that the track parameters are set
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK_EQUAL(ts.predicted(), pc.predicted.parameters());
  BOOST_CHECK(pc.predicted.covariance().has_value());
  BOOST_CHECK_EQUAL(ts.predictedCovariance(), *pc.predicted.covariance());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK_EQUAL(ts.filtered(), pc.filtered.parameters());
  BOOST_CHECK(pc.filtered.covariance().has_value());
  BOOST_CHECK_EQUAL(ts.filteredCovariance(), *pc.filtered.covariance());
  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK_EQUAL(ts.smoothed(), pc.smoothed.parameters());
  BOOST_CHECK(pc.smoothed.covariance().has_value());
  BOOST_CHECK_EQUAL(ts.smoothedCovariance(), *pc.smoothed.covariance());

  // check that the jacobian is set
  BOOST_CHECK(ts.hasJacobian());
  BOOST_CHECK_EQUAL(ts.jacobian(), pc.jacobian);
  BOOST_CHECK_EQUAL(ts.pathLength(), pc.pathLength);
  // check that chi2 is set
  BOOST_CHECK_EQUAL(ts.chi2(), pc.chi2);

  // check that the uncalibrated source link is set
  BOOST_CHECK_EQUAL(&ts.uncalibrated(), &pc.sourceLink);

  // check that the calibrated measurement is set
  BOOST_CHECK(ts.hasCalibrated());
  BOOST_CHECK_EQUAL(&ts.calibratedSourceLink(), &pc.sourceLink);
  BOOST_CHECK_EQUAL(ts.effectiveCalibrated(),
                    pc.sourceLink.parameters.head(nMeasurements));
  BOOST_CHECK_EQUAL(
      ts.effectiveCalibratedCovariance(),
      pc.sourceLink.covariance.topLeftCorner(nMeasurements, nMeasurements));
  {
    ParametersVector mParFull = ParametersVector::Zero();
    CovarianceMatrix mCovFull = CovarianceMatrix::Zero();
    mParFull.head(nMeasurements) = pc.sourceLink.parameters.head(nMeasurements);
    mCovFull.topLeftCorner(nMeasurements, nMeasurements) =
        pc.sourceLink.covariance.topLeftCorner(nMeasurements, nMeasurements);

    auto expMeas = pc.sourceLink.parameters.head(nMeasurements);
    auto expCov =
        pc.sourceLink.covariance.topLeftCorner(nMeasurements, nMeasurements);

    visit_measurement(ts.calibratedSize(), [&](auto N) {
      constexpr size_t measdim = decltype(N)::value;
      BOOST_CHECK_EQUAL(ts.calibrated<measdim>(), expMeas);
      BOOST_CHECK_EQUAL(ts.calibratedCovariance<measdim>(), expCov);
    });
  }

  BOOST_CHECK(ts.hasProjector());
  ActsMatrix<VectorMultiTrajectory::MeasurementSizeMax, eBoundSize> fullProj;
  fullProj.setZero();
  {
    // create a temporary measurement to extract the projector matrix
    auto meas = testSourceLinkCalibratorReturn<VectorMultiTrajectory>(gctx, ts);
    std::visit(
        [&](const auto& m) {
          fullProj.topLeftCorner(nMeasurements, eBoundSize) = m.projector();
        },
        meas);
  }
  BOOST_CHECK_EQUAL(ts.effectiveProjector(),
                    fullProj.topLeftCorner(nMeasurements, eBoundSize));
  BOOST_CHECK_EQUAL(ts.projector(), fullProj);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyAllocations) {
  using namespace Acts::HashedStringLiteral;

  TestTrackState pc(rng, 2u);

  // this should allocate for all components in the trackstate, plus filtered
  VectorMultiTrajectory t;
  size_t i = t.addTrackState(TrackStatePropMask::Predicted |
                             TrackStatePropMask::Filtered |
                             TrackStatePropMask::Jacobian);
  auto tso = t.getTrackState(i);
  fillTrackState(pc, TrackStatePropMask::Predicted, tso);
  fillTrackState(pc, TrackStatePropMask::Filtered, tso);
  fillTrackState(pc, TrackStatePropMask::Jacobian, tso);

  BOOST_CHECK(tso.hasPredicted());
  BOOST_CHECK(tso.hasFiltered());
  BOOST_CHECK(!tso.hasSmoothed());
  BOOST_CHECK(!tso.hasCalibrated());
  BOOST_CHECK(tso.hasJacobian());

  auto tsnone = t.getTrackState(t.addTrackState(TrackStatePropMask::None));
  BOOST_CHECK(!tsnone.has<"predicted"_hash>());
  BOOST_CHECK(!tsnone.has<"filtered"_hash>());
  BOOST_CHECK(!tsnone.has<"smoothed"_hash>());
  BOOST_CHECK(!tsnone.has<"jacobian"_hash>());
  BOOST_CHECK(!tsnone.has<"calibrated"_hash>());
  BOOST_CHECK(!tsnone.has<"projector"_hash>());
  BOOST_CHECK(
      !tsnone.has<"uncalibrated"_hash>());  // separate optional mechanism
  BOOST_CHECK(tsnone.has<"calibratedSourceLink"_hash>());
  BOOST_CHECK(tsnone.has<"referenceSurface"_hash>());
  BOOST_CHECK(tsnone.has<"measdim"_hash>());
  BOOST_CHECK(tsnone.has<"chi2"_hash>());
  BOOST_CHECK(tsnone.has<"pathLength"_hash>());
  BOOST_CHECK(tsnone.has<"typeFlags"_hash>());

  auto tsall = t.getTrackState(t.addTrackState(TrackStatePropMask::All));
  BOOST_CHECK(tsall.has<"predicted"_hash>());
  BOOST_CHECK(tsall.has<"filtered"_hash>());
  BOOST_CHECK(tsall.has<"smoothed"_hash>());
  BOOST_CHECK(tsall.has<"jacobian"_hash>());
  BOOST_CHECK(!tsall.has<"calibrated"_hash>());
  tsall.allocateCalibrated(5);
  BOOST_CHECK(tsall.has<"calibrated"_hash>());
  BOOST_CHECK(tsall.has<"projector"_hash>());
  BOOST_CHECK(!tsall.has<"uncalibrated"_hash>());  // separate optional
                                                   // mechanism: nullptr
  BOOST_CHECK(tsall.has<"calibratedSourceLink"_hash>());
  BOOST_CHECK(tsall.has<"referenceSurface"_hash>());
  BOOST_CHECK(tsall.has<"measdim"_hash>());
  BOOST_CHECK(tsall.has<"chi2"_hash>());
  BOOST_CHECK(tsall.has<"pathLength"_hash>());
  BOOST_CHECK(tsall.has<"typeFlags"_hash>());

  tsall.unset(TrackStatePropMask::Predicted);
  BOOST_CHECK(!tsall.has<"predicted"_hash>());
  tsall.unset(TrackStatePropMask::Filtered);
  BOOST_CHECK(!tsall.has<"filtered"_hash>());
  tsall.unset(TrackStatePropMask::Smoothed);
  BOOST_CHECK(!tsall.has<"smoothed"_hash>());
  tsall.unset(TrackStatePropMask::Jacobian);
  BOOST_CHECK(!tsall.has<"jacobian"_hash>());
  tsall.unset(TrackStatePropMask::Calibrated);
  BOOST_CHECK(!tsall.has<"calibrated"_hash>());
}

BOOST_AUTO_TEST_CASE(TrackStateProxyGetMask) {
  using PM = TrackStatePropMask;

  std::array<PM, 5> values{PM::Predicted, PM::Filtered, PM::Smoothed,
                           PM::Jacobian, PM::Calibrated};
  PM all = std::accumulate(values.begin(), values.end(), PM::None,
                           [](auto a, auto b) { return a | b; });

  VectorMultiTrajectory mj;
  {
    auto ts = mj.getTrackState(mj.addTrackState(PM::All));
    // Calibrated is ignored because we haven't allocated yet
    BOOST_CHECK_EQUAL(ts.getMask(), (all & ~PM::Calibrated));
    ts.allocateCalibrated(4);
    BOOST_CHECK_EQUAL(ts.getMask(), all);
  }
  {
    auto ts = mj.getTrackState(mj.addTrackState(PM::Filtered | PM::Calibrated));
    // Calibrated is ignored because we haven't allocated yet
    BOOST_CHECK_EQUAL(ts.getMask(), PM::Filtered);
    ts.allocateCalibrated(4);
    BOOST_CHECK_EQUAL(ts.getMask(), (PM::Filtered | PM::Calibrated));
  }
  {
    auto ts = mj.getTrackState(
        mj.addTrackState(PM::Filtered | PM::Smoothed | PM::Predicted));
    BOOST_CHECK(ts.getMask() == (PM::Filtered | PM::Smoothed | PM::Predicted));
  }
  {
    for (PM mask : values) {
      auto ts = mj.getTrackState(mj.addTrackState(mask));
      // Calibrated is ignored because we haven't allocated yet
      BOOST_CHECK_EQUAL(ts.getMask(), (mask & ~PM::Calibrated));
    }
  }
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopy) {
  using PM = TrackStatePropMask;

  std::array<PM, 4> values{PM::Predicted, PM::Filtered, PM::Smoothed,
                           PM::Jacobian};

  VectorMultiTrajectory mj;
  auto mkts = [&](PM mask) {
    auto r = mj.getTrackState(mj.addTrackState(mask));
    return r;
  };

  // orthogonal ones
  for (PM a : values) {
    for (PM b : values) {
      auto tsa = mkts(a);
      auto tsb = mkts(b);
      // doesn't work
      if (a != b) {
        BOOST_CHECK_THROW(tsa.copyFrom(tsb), std::runtime_error);
        BOOST_CHECK_THROW(tsb.copyFrom(tsa), std::runtime_error);
      } else {
        tsa.copyFrom(tsb);
        tsb.copyFrom(tsa);
      }
    }
  }

  {
    BOOST_TEST_CHECKPOINT("Calib auto alloc");
    auto tsa = mkts(PM::All);
    auto tsb = mkts(PM::All);
    tsb.allocateCalibrated(5);
    tsb.calibrated<5>().setRandom();
    tsb.calibratedCovariance<5>().setRandom();
    tsa.copyFrom(tsb, PM::All);
    BOOST_CHECK_EQUAL(tsa.calibrated<5>(), tsb.calibrated<5>());
    BOOST_CHECK_EQUAL(tsa.calibratedCovariance<5>(),
                      tsb.calibratedCovariance<5>());
  }

  {
    BOOST_TEST_CHECKPOINT("Copy none");
    auto tsa = mkts(PM::All);
    auto tsb = mkts(PM::All);
    tsa.copyFrom(tsb, PM::None);
  }

  auto ts1 = mkts(PM::Filtered | PM::Predicted);  // this has both
  ts1.filtered().setRandom();
  ts1.filteredCovariance().setRandom();
  ts1.predicted().setRandom();
  ts1.predictedCovariance().setRandom();

  // ((src XOR dst) & src) == 0
  auto ts2 = mkts(PM::Predicted);
  ts2.predicted().setRandom();
  ts2.predictedCovariance().setRandom();

  // they are different before
  BOOST_CHECK(ts1.predicted() != ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() != ts2.predictedCovariance());

  // ts1 -> ts2 fails
  BOOST_CHECK_THROW(ts2.copyFrom(ts1), std::runtime_error);
  BOOST_CHECK(ts1.predicted() != ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() != ts2.predictedCovariance());

  // ts2 -> ts1 is ok
  ts1.copyFrom(ts2);
  BOOST_CHECK(ts1.predicted() == ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() == ts2.predictedCovariance());

  size_t i0 = mj.addTrackState();
  size_t i1 = mj.addTrackState();
  ts1 = mj.getTrackState(i0);
  ts2 = mj.getTrackState(i1);
  TestTrackState rts1(rng, 1u);
  TestTrackState rts2(rng, 2u);
  fillTrackState(rts1, TrackStatePropMask::All, ts1);
  fillTrackState(rts2, TrackStatePropMask::All, ts2);

  auto ots1 = mkts(PM::All);
  auto ots2 = mkts(PM::All);
  // make full copy for later. We prove full copy works right below
  ots1.copyFrom(ts1);
  ots2.copyFrom(ts2);

  BOOST_CHECK_NE(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_NE(ts1.predictedCovariance(), ts2.predictedCovariance());
  BOOST_CHECK_NE(ts1.filtered(), ts2.filtered());
  BOOST_CHECK_NE(ts1.filteredCovariance(), ts2.filteredCovariance());
  BOOST_CHECK_NE(ts1.smoothed(), ts2.smoothed());
  BOOST_CHECK_NE(ts1.smoothedCovariance(), ts2.smoothedCovariance());

  BOOST_CHECK_NE(&ts1.uncalibrated(), &ts2.uncalibrated());

  BOOST_CHECK_NE(&ts1.calibratedSourceLink(), &ts2.calibratedSourceLink());

  visit_measurement(ts1.calibratedSize(), [&](auto N) {
    constexpr size_t measdim = decltype(N)::value;
    BOOST_CHECK_NE(ts1.calibrated<measdim>(), ts2.calibrated<measdim>());
    BOOST_CHECK_NE(ts1.calibratedCovariance<measdim>(),
                   ts2.calibratedCovariance<measdim>());
  });

  BOOST_CHECK_NE(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_NE(ts1.projector(), ts2.projector());

  BOOST_CHECK_NE(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_NE(ts1.chi2(), ts2.chi2());
  BOOST_CHECK_NE(ts1.pathLength(), ts2.pathLength());
  BOOST_CHECK_NE(&ts1.referenceSurface(), &ts2.referenceSurface());

  ts1.copyFrom(ts2);

  BOOST_CHECK_EQUAL(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_EQUAL(ts1.predictedCovariance(), ts2.predictedCovariance());
  BOOST_CHECK_EQUAL(ts1.filtered(), ts2.filtered());
  BOOST_CHECK_EQUAL(ts1.filteredCovariance(), ts2.filteredCovariance());
  BOOST_CHECK_EQUAL(ts1.smoothed(), ts2.smoothed());
  BOOST_CHECK_EQUAL(ts1.smoothedCovariance(), ts2.smoothedCovariance());

  BOOST_CHECK_EQUAL(&ts1.uncalibrated(), &ts2.uncalibrated());

  BOOST_CHECK_EQUAL(&ts1.calibratedSourceLink(), &ts2.calibratedSourceLink());

  visit_measurement(ts1.calibratedSize(), [&](auto N) {
    constexpr size_t measdim = decltype(N)::value;
    BOOST_CHECK_EQUAL(ts1.calibrated<measdim>(), ts2.calibrated<measdim>());
    BOOST_CHECK_EQUAL(ts1.calibratedCovariance<measdim>(),
                      ts2.calibratedCovariance<measdim>());
  });

  BOOST_CHECK_EQUAL(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_EQUAL(ts1.projector(), ts2.projector());

  BOOST_CHECK_EQUAL(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_EQUAL(ts1.chi2(), ts2.chi2());
  BOOST_CHECK_EQUAL(ts1.pathLength(), ts2.pathLength());
  BOOST_CHECK_EQUAL(&ts1.referenceSurface(), &ts2.referenceSurface());

  // full copy proven to work. now let's do partial copy
  ts2 = mkts(PM::Predicted | PM::Jacobian | PM::Calibrated);
  ts2.copyFrom(ots2, PM::Predicted | PM::Jacobian | PM::Calibrated);
  // copy into empty ts, only copy some
  ts1.copyFrom(ots1);  // reset to original
  // is different again
  BOOST_CHECK_NE(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_NE(ts1.predictedCovariance(), ts2.predictedCovariance());

  BOOST_CHECK_NE(&ts1.calibratedSourceLink(), &ts2.calibratedSourceLink());

  visit_measurement(ts1.calibratedSize(), [&](auto N) {
    constexpr size_t measdim = decltype(N)::value;
    BOOST_CHECK_NE(ts1.calibrated<measdim>(), ts2.calibrated<measdim>());
    BOOST_CHECK_NE(ts1.calibratedCovariance<measdim>(),
                   ts2.calibratedCovariance<measdim>());
  });

  BOOST_CHECK_NE(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_NE(ts1.projector(), ts2.projector());

  BOOST_CHECK_NE(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_NE(ts1.chi2(), ts2.chi2());
  BOOST_CHECK_NE(ts1.pathLength(), ts2.pathLength());
  BOOST_CHECK_NE(&ts1.referenceSurface(), &ts2.referenceSurface());

  ts1.copyFrom(ts2);

  // some components are same now
  BOOST_CHECK_EQUAL(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_EQUAL(ts1.predictedCovariance(), ts2.predictedCovariance());

  BOOST_CHECK_EQUAL(&ts1.calibratedSourceLink(), &ts2.calibratedSourceLink());

  visit_measurement(ts1.calibratedSize(), [&](auto N) {
    constexpr size_t measdim = decltype(N)::value;
    BOOST_CHECK_EQUAL(ts1.calibrated<measdim>(), ts2.calibrated<measdim>());
    BOOST_CHECK_EQUAL(ts1.calibratedCovariance<measdim>(),
                      ts2.calibratedCovariance<measdim>());
  });

  BOOST_CHECK_EQUAL(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_EQUAL(ts1.projector(), ts2.projector());

  BOOST_CHECK_EQUAL(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_EQUAL(ts1.chi2(), ts2.chi2());              // always copied
  BOOST_CHECK_EQUAL(ts1.pathLength(), ts2.pathLength());  // always copied
  BOOST_CHECK_EQUAL(&ts1.referenceSurface(),
                    &ts2.referenceSurface());  // always copied
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopyDiffMTJ) {
  using PM = TrackStatePropMask;

  std::array<PM, 4> values{PM::Predicted, PM::Filtered, PM::Smoothed,
                           PM::Jacobian};

  VectorMultiTrajectory mj;
  VectorMultiTrajectory mj2;
  auto mkts = [&](PM mask) {
    auto r = mj.getTrackState(mj.addTrackState(mask));
    return r;
  };
  auto mkts2 = [&](PM mask) {
    auto r = mj2.getTrackState(mj2.addTrackState(mask));
    return r;
  };

  // orthogonal ones
  for (PM a : values) {
    for (PM b : values) {
      auto tsa = mkts(a);
      auto tsb = mkts2(b);
      // doesn't work
      if (a != b) {
        BOOST_CHECK_THROW(tsa.copyFrom(tsb), std::runtime_error);
        BOOST_CHECK_THROW(tsb.copyFrom(tsa), std::runtime_error);
      } else {
        tsa.copyFrom(tsb);
        tsb.copyFrom(tsa);
      }
    }
  }

  // make sure they are actually on different MultiTrajectories
  BOOST_CHECK_EQUAL(mj.size(), values.size() * values.size());
  BOOST_CHECK_EQUAL(mj2.size(), values.size() * values.size());

  auto ts1 = mkts(PM::Filtered | PM::Predicted);  // this has both
  ts1.filtered().setRandom();
  ts1.filteredCovariance().setRandom();
  ts1.predicted().setRandom();
  ts1.predictedCovariance().setRandom();

  // ((src XOR dst) & src) == 0
  auto ts2 = mkts2(PM::Predicted);
  ts2.predicted().setRandom();
  ts2.predictedCovariance().setRandom();

  // they are different before
  BOOST_CHECK(ts1.predicted() != ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() != ts2.predictedCovariance());

  // ts1 -> ts2 fails
  BOOST_CHECK_THROW(ts2.copyFrom(ts1), std::runtime_error);
  BOOST_CHECK(ts1.predicted() != ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() != ts2.predictedCovariance());

  // ts2 -> ts1 is ok
  ts1.copyFrom(ts2);
  BOOST_CHECK(ts1.predicted() == ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() == ts2.predictedCovariance());

  {
    BOOST_TEST_CHECKPOINT("Calib auto alloc");
    auto tsa = mkts(PM::All);
    auto tsb = mkts(PM::All);
    tsb.allocateCalibrated(5);
    tsb.calibrated<5>().setRandom();
    tsb.calibratedCovariance<5>().setRandom();
    tsa.copyFrom(tsb, PM::All);
    BOOST_CHECK_EQUAL(tsa.calibrated<5>(), tsb.calibrated<5>());
    BOOST_CHECK_EQUAL(tsa.calibratedCovariance<5>(),
                      tsb.calibratedCovariance<5>());
  }

  {
    BOOST_TEST_CHECKPOINT("Copy none");
    auto tsa = mkts(PM::All);
    auto tsb = mkts(PM::All);
    tsa.copyFrom(tsb, PM::None);
  }
}

BOOST_AUTO_TEST_CASE(ProxyAssignment) {
  constexpr TrackStatePropMask kMask = TrackStatePropMask::Predicted;
  VectorMultiTrajectory t;
  auto i0 = t.addTrackState(kMask);

  VectorMultiTrajectory::TrackStateProxy tp = t.getTrackState(i0);  // mutable
  VectorMultiTrajectory::TrackStateProxy tp2{tp};       // mutable to mutable
  VectorMultiTrajectory::ConstTrackStateProxy tp3{tp};  // mutable to const
  // const to mutable: this won't compile
  // MultiTrajectory::TrackStateProxy tp4{tp3};
}

// Check if the copy from const does compile, assume the copy is done correctly
BOOST_AUTO_TEST_CASE(CopyFromConst) {
  using PM = TrackStatePropMask;
  VectorMultiTrajectory mj;

  const auto idx_a = mj.addTrackState(PM::All);
  const auto idx_b = mj.addTrackState(PM::All);

  VectorMultiTrajectory::TrackStateProxy mutableProxy = mj.getTrackState(idx_a);

  const VectorMultiTrajectory& cmj = mj;
  VectorMultiTrajectory::ConstTrackStateProxy constProxy =
      cmj.getTrackState(idx_b);

  mutableProxy.copyFrom(constProxy);

  // copy mutable to const: this won't compile
  // constProxy.copyFrom(mutableProxy);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyShare) {
  TestTrackState pc(rng, 2u);

  {
    VectorMultiTrajectory traj;
    size_t ia = traj.addTrackState(TrackStatePropMask::All);
    size_t ib = traj.addTrackState(TrackStatePropMask::None);

    auto tsa = traj.getTrackState(ia);
    auto tsb = traj.getTrackState(ib);

    fillTrackState(pc, TrackStatePropMask::All, tsa);

    BOOST_CHECK(tsa.hasPredicted());
    BOOST_CHECK(!tsb.hasPredicted());
    tsb.shareFrom(tsa, TrackStatePropMask::Predicted);
    BOOST_CHECK(tsa.hasPredicted());
    BOOST_CHECK(tsb.hasPredicted());
    BOOST_CHECK_EQUAL(tsa.predicted(), tsb.predicted());
    BOOST_CHECK_EQUAL(tsa.predictedCovariance(), tsb.predictedCovariance());

    BOOST_CHECK(tsa.hasFiltered());
    BOOST_CHECK(!tsb.hasFiltered());
    tsb.shareFrom(tsa, TrackStatePropMask::Filtered);
    BOOST_CHECK(tsa.hasFiltered());
    BOOST_CHECK(tsb.hasFiltered());
    BOOST_CHECK_EQUAL(tsa.filtered(), tsb.filtered());
    BOOST_CHECK_EQUAL(tsa.filteredCovariance(), tsb.filteredCovariance());

    BOOST_CHECK(tsa.hasSmoothed());
    BOOST_CHECK(!tsb.hasSmoothed());
    tsb.shareFrom(tsa, TrackStatePropMask::Smoothed);
    BOOST_CHECK(tsa.hasSmoothed());
    BOOST_CHECK(tsb.hasSmoothed());
    BOOST_CHECK_EQUAL(tsa.smoothed(), tsb.smoothed());
    BOOST_CHECK_EQUAL(tsa.smoothedCovariance(), tsb.smoothedCovariance());

    BOOST_CHECK(tsa.hasJacobian());
    BOOST_CHECK(!tsb.hasJacobian());
    tsb.shareFrom(tsa, TrackStatePropMask::Jacobian);
    BOOST_CHECK(tsa.hasJacobian());
    BOOST_CHECK(tsb.hasJacobian());
    BOOST_CHECK_EQUAL(tsa.jacobian(), tsb.jacobian());
  }

  {
    VectorMultiTrajectory traj;
    size_t i = traj.addTrackState(TrackStatePropMask::All &
                                  ~TrackStatePropMask::Filtered &
                                  ~TrackStatePropMask::Smoothed);

    auto ts = traj.getTrackState(i);

    BOOST_CHECK(ts.hasPredicted());
    BOOST_CHECK(!ts.hasFiltered());
    BOOST_CHECK(!ts.hasSmoothed());
    ts.predicted().setRandom();
    ts.predictedCovariance().setRandom();

    ts.shareFrom(TrackStatePropMask::Predicted, TrackStatePropMask::Filtered);
    BOOST_CHECK(ts.hasPredicted());
    BOOST_CHECK(ts.hasFiltered());
    BOOST_CHECK(!ts.hasSmoothed());
    BOOST_CHECK_EQUAL(ts.predicted(), ts.filtered());
    BOOST_CHECK_EQUAL(ts.predictedCovariance(), ts.filteredCovariance());

    ts.shareFrom(TrackStatePropMask::Predicted, TrackStatePropMask::Smoothed);
    BOOST_CHECK(ts.hasPredicted());
    BOOST_CHECK(ts.hasFiltered());
    BOOST_CHECK(ts.hasSmoothed());
    BOOST_CHECK_EQUAL(ts.predicted(), ts.filtered());
    BOOST_CHECK_EQUAL(ts.predicted(), ts.smoothed());
    BOOST_CHECK_EQUAL(ts.predictedCovariance(), ts.filteredCovariance());
    BOOST_CHECK_EQUAL(ts.predictedCovariance(), ts.smoothedCovariance());
  }
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryExtraColumns) {
  using namespace HashedStringLiteral;
  using MTJ = VectorMultiTrajectory;

  struct TestColumn {
    double value;
  };

  MTJ traj;
  traj.addColumn<int>("extra_column");
  traj.addColumn<TestColumn>("another_column");

  auto ts1 = traj.getTrackState(traj.addTrackState());
  auto ts2 = traj.getTrackState(
      traj.addTrackState(TrackStatePropMask::All, ts1.index()));
  auto ts3 = traj.getTrackState(
      traj.addTrackState(TrackStatePropMask::All, ts2.index()));

  BOOST_CHECK(ts1.has<"extra_column"_hash>());
  BOOST_CHECK(ts2.has<"extra_column"_hash>());
  BOOST_CHECK(ts3.has<"extra_column"_hash>());

  BOOST_CHECK(ts1.has<"another_column"_hash>());
  BOOST_CHECK(ts2.has<"another_column"_hash>());
  BOOST_CHECK(ts3.has<"another_column"_hash>());

  ts2.component<int, "extra_column"_hash>() = 6;

  BOOST_CHECK_EQUAL((ts2.component<int, "extra_column"_hash>()), 6);

  ts3.component<TestColumn, "another_column"_hash>().value = 7;
  BOOST_CHECK_EQUAL((ts3.component<TestColumn, "another_column"_hash>().value),
                    7);
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryExtraColumnsRuntime) {
  auto runTest = [](auto&& fn) {
    VectorMultiTrajectory mt;
    std::vector<std::string> columns = {"one", "two", "three", "four"};
    for (const auto& c : columns) {
      BOOST_CHECK(!mt.hasColumn(fn(c)));
      mt.addColumn<int>(c);
      BOOST_CHECK(mt.hasColumn(fn(c)));
    }
    for (const auto& c : columns) {
      auto ts1 = mt.getTrackState(mt.addTrackState());
      auto ts2 = mt.getTrackState(mt.addTrackState());
      BOOST_CHECK(ts1.has(fn(c)));
      BOOST_CHECK(ts2.has(fn(c)));
      ts1.component<int>(fn(c)) = 674;
      ts2.component<int>(fn(c)) = 421;
      BOOST_CHECK_EQUAL(ts1.component<int>(fn(c)), 674);
      BOOST_CHECK_EQUAL(ts2.component<int>(fn(c)), 421);
    }
  };

  runTest([](const std::string& c) { return hashString(c.c_str()); });
  // runTest([](const std::string& c) { return c.c_str(); });
  // runTest([](const std::string& c) { return c; });
  // runTest([](std::string_view c) { return c; });
}

BOOST_AUTO_TEST_CASE(MemoryStats) {
  using namespace boost::histogram;
  using cat = axis::category<std::string>;

  VectorMultiTrajectory mt;

  auto stats = mt.statistics();

  std::stringstream ss;
  stats.toStream(ss);
  std::string out = ss.str();
  BOOST_CHECK(!out.empty());
  BOOST_CHECK(out.find("total") != std::string::npos);

  const auto& h = stats.hist;

  auto column_axis = axis::get<cat>(h.axis(0));
  auto type_axis = axis::get<axis::category<>>(h.axis(1));

  for (int t = 0; t < type_axis.size(); t++) {
    for (int c = 0; c < column_axis.size(); c++) {
      BOOST_CHECK_EQUAL(h.at(c, t), 0);
    }
  }

  TestTrackState pc(rng, 2u);
  auto ts = mt.getTrackState(mt.addTrackState());
  fillTrackState(pc, TrackStatePropMask::All, ts);

  stats = mt.statistics();

  for (int t = 0; t < type_axis.size(); t++) {
    BOOST_TEST_CONTEXT((type_axis.bin(t) == 1 ? "meas" : "other"))
    for (int c = 0; c < column_axis.size(); c++) {
      std::string key = column_axis.bin(c);
      BOOST_TEST_CONTEXT("column: " << key) {
        if (t == 0) {
          BOOST_CHECK_NE((h.at(c, t)), 0);
        } else {
          BOOST_CHECK_EQUAL((h.at(c, t)), 0);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
