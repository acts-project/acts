// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE GainMatrixSmoother Tests
#include <boost/optional/optional_io.hpp>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <memory>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace Test {

using Jacobian = BoundParameters::CovMatrix_t;
using Covariance = BoundSymMatrix;

using SourceLink = MinimalSourceLink;

template <ParID_t... params>
using MeasurementType = Measurement<SourceLink, params...>;
using TrackState = TrackState<SourceLink, BoundParameters>;

// Create a test context
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_CASE(gain_matrix_smoother) {
  // Make dummy measurement
  auto plane1 = Surface::makeShared<PlaneSurface>(Vector3D::UnitX() * 1,
                                                  Vector3D::UnitX());
  auto plane2 = Surface::makeShared<PlaneSurface>(Vector3D::UnitX() * 2,
                                                  Vector3D::UnitX());
  auto plane3 = Surface::makeShared<PlaneSurface>(Vector3D::UnitX() * 3,
                                                  Vector3D::UnitX());

  ActsSymMatrixD<2> cov;
  cov << 0.04, 0, 0, 0.1;
  FittableMeasurement<SourceLink> meas1(
      MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1>(
          plane1, {}, std::move(cov), -0.1, 0.45));

  cov << 0.04, 0, 0, 0.1;
  FittableMeasurement<SourceLink> meas2(
      MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1>(
          plane2, {}, std::move(cov), -0.2, 0.35));

  cov << 0.04, 0, 0, 0.1;
  FittableMeasurement<SourceLink> meas3(
      MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1>(
          plane3, {}, std::move(cov), -0.05, 0.25));

  std::vector<TrackState> trackStates;

  // Make dummy track parameter
  Covariance covTrk;
  covTrk.setIdentity();
  covTrk.diagonal() << 0.08, 0.3, 1, 1, 1, 1;
  // covTrk << 0.08, 0, 0, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
  // 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
  BoundVector parValues;
  parValues << 0.3, 0.5, 0.5 * M_PI, 0., 1 / 100., 0.;
  BoundParameters pars(
      tgContext,
      std::make_unique<const BoundParameters::CovMatrix_t>(std::move(covTrk)),
      parValues, plane1);

  parValues << 0.301, 0.503, 0.5 * M_PI, 0., 1 / 100., 0.;
  BoundParameters filt(
      tgContext,
      std::make_unique<const BoundParameters::CovMatrix_t>(std::move(covTrk)),
      parValues, plane1);

  TrackState ts{SourceLink{&meas1}};

  ts.parameter.predicted = std::move(pars);
  ts.parameter.filtered = std::move(filt);
  ts.parameter.jacobian = Jacobian::Identity();
  ts.parameter.pathLength = 1.;

  trackStates.push_back(ts);

  ts = TrackState{SourceLink{&meas2}};

  parValues << 0.2, 0.5, 0.5 * M_PI, 0., 1 / 100., 0.;
  pars = BoundParameters(
      tgContext,
      std::make_unique<const BoundParameters::CovMatrix_t>(std::move(covTrk)),
      parValues, plane2);
  parValues << 0.27, 0.53, 0.5 * M_PI, 0., 1 / 100., 0.;
  filt = BoundParameters(
      tgContext,
      std::make_unique<const BoundParameters::CovMatrix_t>(std::move(covTrk)),
      parValues, plane2);

  ts.parameter.predicted = std::move(pars);
  ts.parameter.filtered = std::move(filt);
  ts.parameter.jacobian = Jacobian::Identity();
  ts.parameter.pathLength = 2.;

  trackStates.push_back(ts);

  ts = TrackState{SourceLink{&meas3}};

  parValues << 0.35, 0.49, 0.5 * M_PI, 0., 1 / 100., 0.;
  pars = BoundParameters(
      tgContext,
      std::make_unique<const BoundParameters::CovMatrix_t>(std::move(covTrk)),
      parValues, plane3);
  parValues << 0.33, 0.43, 0.5 * M_PI, 0., 1 / 100., 0.;
  filt = BoundParameters(
      tgContext,
      std::make_unique<const BoundParameters::CovMatrix_t>(std::move(covTrk)),
      parValues, plane3);

  ts.parameter.predicted = std::move(pars);
  ts.parameter.filtered = std::move(filt);
  ts.parameter.jacobian = Jacobian::Identity();
  ts.parameter.pathLength = 3.;

  trackStates.push_back(ts);

  // "smooth" these three track states

  GainMatrixSmoother<BoundParameters> gms;
  BOOST_CHECK(gms(tgContext, trackStates));

  // Regression tests, only tests very basic correctness of the math, but tests
  // for regressions in the result.

  auto& ts1 = trackStates.at(0);
  BOOST_CHECK(ts1.parameter.smoothed);
  BOOST_CHECK_NE(ts1.parameter.filtered->parameters(),
                 ts1.parameter.smoothed->parameters());

  double tol = 1e-6;

  BoundVector expPars;
  expPars << 0.3510000, 0.4730000, 1.5707963, 0.0000000, 0.0100000, 0.0000000;
  CHECK_CLOSE_ABS(ts1.parameter.smoothed->parameters(), expPars, tol);
  Covariance expCov;
  expCov.setIdentity();
  expCov.diagonal() << 0.0800000, 0.3000000, 1.0000000, 1.0000000, 1.0000000,
      1.0000000;
  CHECK_CLOSE_ABS(*ts1.parameter.smoothed->covariance(), expCov, tol);

  auto& ts2 = trackStates.at(1);
  BOOST_CHECK(ts2.parameter.smoothed);
  BOOST_CHECK_NE(ts2.parameter.filtered->parameters(),
                 ts2.parameter.smoothed->parameters());

  expPars << 0.2500000, 0.4700000, 1.5707963, 0.0000000, 0.0100000, 0.0000000;
  CHECK_CLOSE_ABS(ts2.parameter.smoothed->parameters(), expPars, tol);
  CHECK_CLOSE_ABS(*ts2.parameter.smoothed->covariance(), expCov, tol);

  auto& ts3 = trackStates.at(2);
  BOOST_CHECK(ts3.parameter.smoothed);
  // last one, smoothed == filtered
  BOOST_CHECK_EQUAL(ts3.parameter.filtered->parameters(),
                    ts3.parameter.smoothed->parameters());

  expPars << 0.3300000, 0.4300000, 1.5707963, 0.0000000, 0.0100000, 0.0000000;
  CHECK_CLOSE_ABS(ts3.parameter.smoothed->parameters(), expPars, tol);
  CHECK_CLOSE_ABS(*ts3.parameter.smoothed->covariance(), expCov, tol);
}

}  // namespace Test
}  // namespace Acts
