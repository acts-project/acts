// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/EDM4hep/EDM4hepUtil.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;
BOOST_AUTO_TEST_SUITE(EDM4hepParameterConversion)

BOOST_AUTO_TEST_CASE(ConvertTrackParametersToEdm4hepWithPerigee) {
  auto refSurface = Surface::makeShared<PerigeeSurface>(Vector3{50, 30, 20});

  BoundVector par;
  par << 1_mm, 5_mm, 0, M_PI_2, -1 / 1_GeV,
      5_ns;  // -> perpendicular to perigee and pointing right, should be PCA

  BoundMatrix cov;
  cov.setIdentity();

  SingleBoundTrackParameters<SinglyCharged> boundPar{refSurface, par, cov};

  double Bz = 2_T;

  Acts::GeometryContext gctx;

  EDM4hepUtil::detail::Parameters converted =
      EDM4hepUtil::detail::convertTrackParametersToEdm4hep(gctx, Bz, boundPar);

  BOOST_CHECK(converted.covariance.has_value());
  BOOST_CHECK(converted.surface);

  // input is already on perigee, should not be modified
  BOOST_CHECK_EQUAL(par.template head<2>(),
                    converted.values.template head<2>());
  BOOST_CHECK_EQUAL(
      (converted.covariance.value().template topLeftCorner<4, 4>()),
      ActsSymMatrix<4>::Identity());
  BOOST_CHECK(converted.covariance.value()(4, 4) > 0);

  // convert back for roundtrip test

  SingleBoundTrackParameters<SinglyCharged> roundtripPar =
      EDM4hepUtil::detail::convertTrackParametersFromEdm4hep<SinglyCharged>(
          Bz, converted);

  BOOST_CHECK(roundtripPar.parameters().isApprox(boundPar.parameters()));
  BOOST_CHECK(roundtripPar.covariance().value().isApprox(
      boundPar.covariance().value()));
}

BOOST_AUTO_TEST_CASE(ConvertTrackParametersToEdm4hepWithOutPerigee) {
  auto refSurface = Surface::makeShared<PlaneSurface>(
      Vector3{50, 30, 20}, Vector3{1, 1, 0.3}.normalized());

  BoundVector par;
  par << 1_mm, 5_mm, M_PI / 4., M_PI_2, -1 / 1_GeV, 5_ns;

  BoundMatrix cov;
  cov.setIdentity();

  SingleBoundTrackParameters<SinglyCharged> boundPar{refSurface, par, cov};

  double Bz = 2_T;

  Acts::GeometryContext gctx;

  EDM4hepUtil::detail::Parameters converted =
      EDM4hepUtil::detail::convertTrackParametersToEdm4hep(gctx, Bz, boundPar);

  BOOST_CHECK(converted.covariance.has_value());
  BOOST_CHECK(converted.surface);

  // input is not a perigee, so new params should be at 0, 0 on ad-hoc perigee
  BOOST_CHECK_EQUAL(converted.values.template head<2>(), (Vector2{0, 0}));
  BOOST_CHECK_EQUAL(converted.values[2], par[2]);

  // std::cout << converted.covariance.value() << std::endl;

  BOOST_CHECK_EQUAL(
      (converted.covariance.value().template topLeftCorner<4, 4>()),
      ActsSymMatrix<4>::Identity());
  BOOST_CHECK(converted.covariance.value()(4, 4) > 0);

  // convert back for roundtrip test
  SingleBoundTrackParameters<SinglyCharged> roundtripPar =
      EDM4hepUtil::detail::convertTrackParametersFromEdm4hep<SinglyCharged>(
          Bz, converted);

  BOOST_CHECK_EQUAL(roundtripPar.parameters().template head<2>(),
                    (Vector2{0, 0}));
  BOOST_CHECK_EQUAL(roundtripPar.parameters().tail<4>(), par.tail<4>());
  BOOST_CHECK(roundtripPar.covariance().value().isApprox(
      boundPar.covariance().value()));
}

BOOST_AUTO_TEST_SUITE_END()
