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
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackHelpers.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/EDM4hep/EDM4hepUtil.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <algorithm>
#include <random>

#include <edm4hep/TrackCollection.h>

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
  cov(5, 5) = 25_ns;

  BoundTrackParameters boundPar{refSurface, par, cov,
                                ParticleHypothesis::pion()};

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
      ActsSquareMatrix<4>::Identity());
  BOOST_CHECK_GT(converted.covariance.value()(4, 4), 0);
  BOOST_CHECK_EQUAL(converted.covariance.value()(5, 5), 25_ns);

  // convert back for roundtrip test

  BoundTrackParameters roundtripPar =
      EDM4hepUtil::detail::convertTrackParametersFromEdm4hep(Bz, converted);

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
  cov(5, 5) = 25_ns;

  BoundTrackParameters boundPar{refSurface, par, cov,
                                ParticleHypothesis::pion()};

  double Bz = 2_T;

  Acts::GeometryContext gctx;

  EDM4hepUtil::detail::Parameters converted =
      EDM4hepUtil::detail::convertTrackParametersToEdm4hep(gctx, Bz, boundPar);

  BOOST_CHECK(converted.covariance.has_value());
  BOOST_CHECK(converted.surface);

  // input is not a perigee, so new params should be at 0, 0 on ad-hoc perigee
  BOOST_CHECK_EQUAL(converted.values.template head<2>(), (Vector2{0, 0}));
  CHECK_CLOSE_ABS(converted.values[2], par[2], 1e-6);

  BOOST_CHECK((converted.covariance.value().template topLeftCorner<4, 4>())
                  .isApprox(ActsSquareMatrix<4>::Identity()));
  BOOST_CHECK_GT(converted.covariance.value()(4, 4), 0);
  BOOST_CHECK_EQUAL(converted.covariance.value()(5, 5), 25_ns);

  // convert back for roundtrip test
  BoundTrackParameters roundtripPar =
      EDM4hepUtil::detail::convertTrackParametersFromEdm4hep(Bz, converted);

  BOOST_CHECK_EQUAL(roundtripPar.parameters().template head<2>(),
                    (Vector2{0, 0}));
  BOOST_CHECK(roundtripPar.parameters().tail<4>().isApprox(par.tail<4>()));
  BOOST_CHECK(roundtripPar.covariance().value().isApprox(
      boundPar.covariance().value()));
}

BOOST_AUTO_TEST_CASE(ConvertTrackParametersToEdm4hepWithPerigeeNoCov) {
  auto refSurface = Surface::makeShared<PerigeeSurface>(Vector3{50, 30, 20});

  BoundVector par;
  par << 1_mm, 5_mm, 0, M_PI_2, -1 / 1_GeV,
      5_ns;  // -> perpendicular to perigee and pointing right, should be PCA

  BoundTrackParameters boundPar{refSurface, par, std::nullopt,
                                ParticleHypothesis::pion()};

  double Bz = 2_T;

  Acts::GeometryContext gctx;

  EDM4hepUtil::detail::Parameters converted =
      EDM4hepUtil::detail::convertTrackParametersToEdm4hep(gctx, Bz, boundPar);

  BOOST_CHECK(!converted.covariance.has_value());
  BOOST_CHECK(converted.surface);

  // input is already on perigee, should not be modified
  BOOST_CHECK_EQUAL(par.template head<2>(),
                    converted.values.template head<2>());

  // convert back for roundtrip test

  BoundTrackParameters roundtripPar =
      EDM4hepUtil::detail::convertTrackParametersFromEdm4hep(Bz, converted);

  BOOST_CHECK(roundtripPar.parameters().isApprox(boundPar.parameters()));
  BOOST_CHECK(!roundtripPar.covariance().has_value());
}

BOOST_AUTO_TEST_CASE(ConvertTrackParametersToEdm4hepWithOutPerigeeNoCov) {
  auto refSurface = Surface::makeShared<PlaneSurface>(
      Vector3{50, 30, 20}, Vector3{1, 1, 0.3}.normalized());

  BoundVector par;
  par << 1_mm, 5_mm, M_PI / 4., M_PI_2, -1 / 1_GeV, 5_ns;

  BoundTrackParameters boundPar{refSurface, par, std::nullopt,
                                ParticleHypothesis::pion()};

  double Bz = 2_T;

  Acts::GeometryContext gctx;

  EDM4hepUtil::detail::Parameters converted =
      EDM4hepUtil::detail::convertTrackParametersToEdm4hep(gctx, Bz, boundPar);

  BOOST_CHECK(!converted.covariance.has_value());
  BOOST_CHECK(converted.surface);

  // input is not a perigee, so new params should be at 0, 0 on ad-hoc perigee
  BOOST_CHECK_EQUAL(converted.values.template head<2>(), (Vector2{0, 0}));
  CHECK_CLOSE_ABS(converted.values[2], par[2], 1e-6);

  // convert back for roundtrip test
  BoundTrackParameters roundtripPar =
      EDM4hepUtil::detail::convertTrackParametersFromEdm4hep(Bz, converted);

  BOOST_CHECK_EQUAL(roundtripPar.parameters().template head<2>(),
                    (Vector2{0, 0}));
  BOOST_CHECK(roundtripPar.parameters().tail<4>().isApprox(par.tail<4>()));
  BOOST_CHECK(!roundtripPar.covariance().has_value());
}

BOOST_AUTO_TEST_CASE(CovariancePacking) {
  BoundMatrix m;
  // clang-format off
  m << 1, 2, 3, 4, 5, 6,
       2, 2, 3, 4, 5, 6,
       3, 3, 3, 4, 5, 6,
       4, 4, 4, 4, 5, 6,
       5, 5, 5, 5, 5, 6,
       6, 6, 6, 6, 6, 6;
  // clang-format on

  std::array<float, 21> values{};
  EDM4hepUtil::detail::packCovariance(m, values.data());

  BoundMatrix m2;
  m2.setZero();
  EDM4hepUtil::detail::unpackCovariance(values.data(), m2);

  CHECK_CLOSE_ABS(m, m2, 1e-9);
}

BOOST_AUTO_TEST_CASE(RoundTripTests) {
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  using mutable_proxy_t = decltype(tracks)::TrackProxy;
  using const_proxy_t = decltype(tracks)::ConstTrackProxy;

  std::mt19937 rng{42};
  std::normal_distribution<double> gauss(0., 1.);
  std::uniform_real_distribution<double> f(-1, 1);
  std::uniform_real_distribution<double> r(0, 1);
  std::uniform_int_distribution<unsigned int> nTracks(2, 20);
  std::uniform_int_distribution<unsigned int> nTs(1, 20);
  std::uniform_real_distribution<double> phiDist(-M_PI, M_PI);
  std::uniform_real_distribution<double> etaDist(-4, 4);
  std::uniform_real_distribution<double> ptDist(1_MeV, 10_GeV);
  std::uniform_real_distribution<double> qDist(0., 1.);

  auto genParams = [&]() -> std::pair<BoundVector, BoundMatrix> {
    double d0 = 20_um * gauss(rng);
    double z0 = 20_mm * gauss(rng);
    double phi = phiDist(rng);
    double eta = etaDist(rng);
    double theta = 2 * atan(exp(-eta));
    double pt = ptDist(rng);
    double p = pt / sin(theta);
    double charge = qDist(rng) > 0.5 ? 1. : -1.;
    double qop = charge / p;
    double t = 5_ns * gauss(rng);

    BoundVector par;
    par << d0, z0, phi, theta, qop, t;
    BoundMatrix cov;
    cov = BoundMatrix::Identity();
    cov.diagonal() << 20_um * 20_um, 20_mm * 20_mm, 0.1, 0.1, 1_GeV, 25_ns;
    return {par, cov};
  };

  size_t numT = nTracks(rng);
  for (size_t t = 0; t < numT; t++) {
    auto track = tracks.getTrack(tracks.addTrack());
    {
      auto [par, cov] = genParams();
      track.parameters() = par;
      track.covariance() = cov;
    }
    track.setReferenceSurface(
        Acts::Surface::makeShared<PerigeeSurface>(Vector3{0, 0, 0}));

    size_t numTs = nTs(rng);
    for (size_t i = 0; i < numTs; i++) {
      auto ts = track.appendTrackState(TrackStatePropMask::Smoothed);
      double crit = r(rng);
      if (crit < 0.1) {
        ts.typeFlags().set(TrackStateFlag::HoleFlag);
        continue;
      } else if (crit < 0.2) {
        ts.typeFlags().set(TrackStateFlag::OutlierFlag);
        continue;
      } else if (crit < 0.3) {
        ts.typeFlags().set(TrackStateFlag::SharedHitFlag);
      } else if (crit < 0.4) {
        ts.typeFlags().set(TrackStateFlag::MaterialFlag);
        continue;
      }

      ts.typeFlags().set(TrackStateFlag::MeasurementFlag);

      auto [par, cov] = genParams();
      ts.smoothed() = par;
      ts.smoothedCovariance() = cov;
      Vector3 pos;
      pos << 1000 * f(rng), 1000 * f(rng), 3000 * f(rng);
      ts.setReferenceSurface(Acts::Surface::makeShared<PerigeeSurface>(pos));
    }

    calculateTrackQuantities(track);
  }

  edm4hep::TrackCollection edm4hepTracks;

  Acts::GeometryContext gctx;

  double Bz = 3_T;

  auto logger = getDefaultLogger("EDM4hep", Logging::INFO);

  for (const_proxy_t track : tracks) {
    auto to = edm4hepTracks.create();
    EDM4hepUtil::writeTrack(gctx, track, to, Bz, *logger);
  }

  BOOST_CHECK_EQUAL(edm4hepTracks.size(), tracks.size());

  auto tIt = tracks.begin();
  for (auto edm4hepTrack : edm4hepTracks) {
    auto track = *tIt;
    BOOST_CHECK_EQUAL(track.nMeasurements(),
                      edm4hepTrack.trackStates_size() - 1);

    ++tIt;
  }

  const edm4hep::TrackCollection& edm4hepTracksConst = edm4hepTracks;

  TrackContainer readTracks(std::make_shared<Acts::VectorTrackContainer>(),
                            std::make_shared<Acts::VectorMultiTrajectory>());

  for (const auto edm4hepTrack : edm4hepTracksConst) {
    EDM4hepUtil::readTrack(
        edm4hepTrack, readTracks.getTrack(readTracks.addTrack()), Bz, *logger);
  }

  BOOST_CHECK_EQUAL(tracks.size(), readTracks.size());
  size_t t = 0;

  auto origTrackIt = tracks.begin();
  auto readTrackIt = readTracks.begin();
  while (origTrackIt != tracks.end() && readTrackIt != readTracks.end()) {
    BOOST_TEST_INFO_SCOPE("Track #" << t);
    auto orig = *origTrackIt;
    auto read = *readTrackIt;

    CHECK_CLOSE_OR_SMALL(orig.parameters(), read.parameters(), 1e-5, 1e-8);
    CHECK_CLOSE_OR_SMALL(orig.covariance(), read.covariance(), 1e-5, 1e-8);
    BOOST_CHECK_EQUAL(orig.referenceSurface().center(gctx),
                      read.referenceSurface().center(gctx));

    auto origTsIt = orig.trackStatesReversed().begin();
    auto readTsIt = read.trackStatesReversed().begin();

    size_t tsi = 0;
    while (origTsIt != orig.trackStatesReversed().end() &&
           readTsIt != read.trackStatesReversed().end()) {
      BOOST_TEST_INFO_SCOPE("TS: #" << tsi);
      auto nextMeas = std::find_if(
          origTsIt, orig.trackStatesReversed().end(), [](const auto& ts) {
            return ts.typeFlags().test(TrackStateFlag::MeasurementFlag);
          });
      BOOST_CHECK(nextMeas != orig.trackStatesReversed().end());
      origTsIt = nextMeas;
      auto origTs = *origTsIt;
      auto readTs = *readTsIt;

      BOOST_TEST_INFO_SCOPE(
          "orig parameters: " << origTs.parameters().transpose());
      BOOST_TEST_INFO_SCOPE(
          "read parameters: " << readTs.parameters().transpose());
      CHECK_CLOSE_OR_SMALL(origTs.smoothedCovariance(),
                           readTs.smoothedCovariance(), 1e-5, 1e-6);
      Vector3 newCenter = readTs.referenceSurface().center(
          gctx);  // new center is a perigee, but should be on the other
                  // surface
      BOOST_CHECK(origTs.referenceSurface().isOnSurface(gctx, newCenter,
                                                        Vector3::Zero()));

      // global hit positions should be the same
      Vector3 readGlobal = readTs.referenceSurface().localToGlobal(
          gctx, readTs.parameters().template head<2>(), Vector3::Zero());
      Vector3 origGlobal = origTs.referenceSurface().localToGlobal(
          gctx, origTs.parameters().template head<2>(), Vector3::Zero());
      CHECK_CLOSE_ABS(readGlobal, origGlobal, 1e-3);
      ++origTsIt;
      ++readTsIt;
      tsi++;
    }
    ++origTrackIt;
    ++readTrackIt;

    t++;
  }
}

BOOST_AUTO_TEST_SUITE_END()
