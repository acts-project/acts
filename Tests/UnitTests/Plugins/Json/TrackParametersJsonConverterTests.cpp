// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsPlugins/Json/TrackParametersJsonConverter.hpp"

#include <memory>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(TrackParametersJsonIO) {
  GeometryContext gctx;

  // Track parameters
  Vector4 position(1., 2., 3., 4.);
  double phi = 0.1;
  double theta = 0.2;
  double qOverP = 3.0;
  ParticleHypothesis particle = ParticleHypothesis::electron();
  FreeMatrix freeCov = FreeMatrix::Identity();
  BoundMatrix boundCov = BoundMatrix::Identity();

  auto surface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(10., 10.));
  surface->assignGeometryId(GeometryIdentifier(1u));

  // Free track parameters conversion
  FreeTrackParameters ftp(position, phi, theta, qOverP, freeCov, particle);

  nlohmann::json ftpJson = ftp;

  FreeTrackParameters ftpRead = ftpJson;

  BOOST_CHECK_EQUAL(ftp.position(), ftpRead.position());
  BOOST_CHECK_EQUAL(ftp.direction(), ftpRead.direction());
  BOOST_CHECK_EQUAL(ftp.qOverP(), ftpRead.qOverP());
  BOOST_CHECK_EQUAL(ftp.covariance().value(), ftpRead.covariance().value());
  BOOST_CHECK_EQUAL(ftp.particleHypothesis(), ftpRead.particleHypothesis());

  // Curvilinear track parameters conversion
  BoundTrackParameters ctp = BoundTrackParameters::createCurvilinear(
      position, phi, theta, qOverP, boundCov, particle);

  nlohmann::json ctpJson = ctp;

  BoundTrackParameters ctpRead = ctpJson;

  BOOST_CHECK_EQUAL(ctp.position(gctx), ctpRead.position(gctx));
  BOOST_CHECK_EQUAL(ctp.direction(), ctpRead.direction());
  BOOST_CHECK_EQUAL(ctp.qOverP(), ctpRead.qOverP());
  BOOST_CHECK_EQUAL(ctp.covariance().value(), ctpRead.covariance().value());
  BOOST_CHECK_EQUAL(ctp.particleHypothesis(), ctpRead.particleHypothesis());

  BOOST_CHECK(ctp.referenceSurface().localToGlobal(gctx).isApprox(
      ctpRead.referenceSurface().localToGlobal(gctx)));
  BOOST_CHECK_EQUAL(ctp.referenceSurface().geometryId(),
                    ctpRead.referenceSurface().geometryId());
  BOOST_CHECK_EQUAL(ctp.referenceSurface().bounds(),
                    ctpRead.referenceSurface().bounds());

  // Bound track parameters conversion
  BoundVector boundPosition{1., 2., 3., 4., 5., 6.};
  BoundTrackParameters btp(surface, boundPosition, boundCov, particle);

  nlohmann::json btpJson = btp;

  BoundTrackParameters btpRead = btpJson;

  BOOST_CHECK_EQUAL(btp.position(gctx), btpRead.position(gctx));
  BOOST_CHECK_EQUAL(btp.direction(), btpRead.direction());
  BOOST_CHECK_EQUAL(btp.qOverP(), btpRead.qOverP());
  BOOST_CHECK_EQUAL(btp.covariance().value(), btpRead.covariance().value());
  BOOST_CHECK_EQUAL(btp.particleHypothesis(), btpRead.particleHypothesis());

  BOOST_CHECK(btp.referenceSurface().localToGlobal(gctx).isApprox(
      btpRead.referenceSurface().localToGlobal(gctx)));
  BOOST_CHECK_EQUAL(btp.referenceSurface().geometryId(),
                    btpRead.referenceSurface().geometryId());
  BOOST_CHECK_EQUAL(btp.referenceSurface().bounds(),
                    btpRead.referenceSurface().bounds());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
