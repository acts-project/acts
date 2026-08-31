// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <array>
#include <numbers>
#include <random>

#include <edm4hep/Constants.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackState.h>
#include <edm4hep/TrackerHit3DCollection.h>

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EDM4HepSuite)

BOOST_AUTO_TEST_CASE(JacobianRoundtrip) {
  BoundVector par;
  par << 1_mm, 5_mm, 0.1, std::numbers::pi / 2. * 0.9, -1 / 1_GeV, 5_ns;

  BoundMatrix cov;
  cov.setIdentity();

  double Bz = 2_T;

  double tanLambda = std::tan(std::numbers::pi / 2. - par[eBoundTheta]);
  double omega = par[eBoundQOverP] / std::sin(par[eBoundTheta]) * Bz;

  auto J1 = EDM4hepUtil::detail::jacobianToEdm4hep(par[eBoundTheta],
                                                   par[eBoundQOverP], Bz);

  BoundMatrix cov2 = J1 * cov * J1.transpose();

  auto J2 = EDM4hepUtil::detail::jacobianFromEdm4hep(tanLambda, omega, Bz);

  BoundMatrix cov3 = J2 * cov2 * J2.transpose();

  CHECK_CLOSE_ABS(cov, cov3, 1e-9);
}

BOOST_AUTO_TEST_CASE(ConvertTrackParametersToEdm4hepWithPerigee) {
  auto refSurface = Surface::makeShared<PerigeeSurface>(Vector3{50, 30, 20});

  BoundVector par;
  par << 1_mm, 5_mm, 0, std::numbers::pi / 2., -1 / 1_GeV,
      5_ns;  // -> perpendicular to perigee and pointing right, should be PCA

  BoundMatrix cov;
  cov.setIdentity();
  cov(5, 5) = 25_ns;

  BoundTrackParameters boundPar{refSurface, par, cov,
                                ParticleHypothesis::pion()};

  double Bz = 2_T;

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  EDM4hepUtil::detail::Parameters converted =
      EDM4hepUtil::detail::convertTrackParametersToEdm4hep(gctx, Bz, boundPar);

  BOOST_CHECK(converted.covariance.has_value());
  BOOST_CHECK(converted.surface);

  // input is already on perigee, should not be modified
  BOOST_CHECK_EQUAL(par.template head<2>(),
                    converted.values.template head<2>());
  BOOST_CHECK_EQUAL(
      (converted.covariance.value().template topLeftCorner<4, 4>()),
      SquareMatrix<4>::Identity());
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
  std::shared_ptr<PlaneSurface> planeSurface =
      CurvilinearSurface(Vector3{50, 30, 20}, Vector3{1, 1, 0.3}.normalized())
          .planeSurface();

  BoundVector par;
  par << 1_mm, 5_mm, std::numbers::pi / 4., std::numbers::pi / 2. * 0.9,
      -1 / 1_GeV, 5_ns;

  BoundMatrix cov;
  cov.setIdentity();
  cov(5, 5) = 25_ns;

  BoundTrackParameters planePar{planeSurface, par, cov,
                                ParticleHypothesis::pion()};

  double Bz = 2_T;

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  EDM4hepUtil::detail::Parameters converted =
      EDM4hepUtil::detail::convertTrackParametersToEdm4hep(gctx, Bz, planePar);

  BOOST_CHECK(converted.covariance.has_value());
  BOOST_CHECK(converted.surface);

  // input is not a perigee, so new params should be at 0, 0 on ad-hoc perigee
  BOOST_CHECK_EQUAL(converted.values.template head<2>(), (Vector2{0, 0}));
  CHECK_CLOSE_ABS(converted.values[2], par[2], 1e-6);

  BOOST_CHECK_EQUAL(converted.covariance.value()(0, 0), 1);

  BOOST_CHECK_LT(converted.covariance.value()(1, 1), 1.2);
  BOOST_CHECK_GT(converted.covariance.value()(1, 1), 1);

  CHECK_CLOSE_ABS(converted.covariance.value()(2, 2), 1, 1e-6);

  BOOST_CHECK_GT(converted.covariance.value()(3, 3), 1);
  BOOST_CHECK_LT(converted.covariance.value()(3, 3), 1.2);

  BOOST_CHECK_GT(converted.covariance.value()(4, 4), 0);
  BOOST_CHECK_EQUAL(converted.covariance.value()(5, 5), 25_ns);

  // convert back for roundtrip test
  BoundTrackParameters roundtripPar =
      EDM4hepUtil::detail::convertTrackParametersFromEdm4hep(Bz, converted);

  BOOST_CHECK_NE(
      dynamic_cast<const PerigeeSurface*>(&roundtripPar.referenceSurface()),
      nullptr);

  BOOST_CHECK((converted.covariance.value().topLeftCorner<3, 3>().isApprox(
      roundtripPar.covariance().value().topLeftCorner<3, 3>())));
  CHECK_CLOSE_ABS(roundtripPar.covariance().value()(3, 3), 1, 1e-6);
  CHECK_CLOSE_ABS(roundtripPar.covariance().value()(4, 4), 1, 1e-6);
  BOOST_CHECK_EQUAL(roundtripPar.covariance().value()(5, 5), 25_ns);

  auto roundtripPlaneBoundParams =
      detail::boundToBoundConversion(gctx, roundtripPar, *planeSurface).value();

  BOOST_CHECK(roundtripPlaneBoundParams.parameters().isApprox(par));

  CHECK_CLOSE_COVARIANCE(roundtripPlaneBoundParams.covariance().value(),
                         planePar.covariance().value(), 1e-3);
}

BOOST_AUTO_TEST_CASE(ConvertTrackParametersToEdm4hepWithPerigeeNoCov) {
  auto refSurface = Surface::makeShared<PerigeeSurface>(Vector3{50, 30, 20});

  BoundVector par;
  par << 1_mm, 5_mm, 0, std::numbers::pi / 2., -1 / 1_GeV,
      5_ns;  // -> perpendicular to perigee and pointing right, should be PCA

  BoundTrackParameters boundPar{refSurface, par, std::nullopt,
                                ParticleHypothesis::pion()};

  double Bz = 2_T;

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

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
  std::shared_ptr<PlaneSurface> refSurface =
      CurvilinearSurface(Vector3{50, 30, 20}, Vector3{1, 1, 0.3}.normalized())
          .planeSurface();

  BoundVector par;
  par << 1_mm, 5_mm, std::numbers::pi / 4., std::numbers::pi / 2., -1 / 1_GeV,
      5_ns;

  BoundTrackParameters boundPar{refSurface, par, std::nullopt,
                                ParticleHypothesis::pion()};

  double Bz = 2_T;

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

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

  edm4hep::CovMatrix6f values{};
  EDM4hepUtil::detail::packCovariance(m, values);

  BoundMatrix m2;
  m2.setZero();
  EDM4hepUtil::detail::unpackCovariance(values, m2);

  CHECK_CLOSE_ABS(m, m2, 1e-9);
}

// The packed covariance has to follow the EDM4hep parameter order
// (edm4hep::TrackParams: d0, phi, omega, z0, tanLambda, time), which is a
// different permutation from the order used internally by
// detail::Parameters::values (d0, z0, phi, tanLambda, omega, time). A
// pack/unpack round trip cannot catch a mismatch here, since it would be
// symmetric, so check the absolute slots against EDM4hep's own accessor.
BOOST_AUTO_TEST_CASE(CovariancePackingConvention) {
  // values index -> edm4hep parameter it represents
  constexpr std::array<edm4hep::TrackParams, 6> params{
      edm4hep::TrackParams::d0,         // values[0]
      edm4hep::TrackParams::z0,         // values[1]
      edm4hep::TrackParams::phi,        // values[2]
      edm4hep::TrackParams::tanLambda,  // values[3]
      edm4hep::TrackParams::omega,      // values[4]
      edm4hep::TrackParams::time,       // values[5]
  };

  // Distinct, symmetric entries so every slot is identifiable
  SquareMatrix<6> m;
  for (std::size_t i = 0; i < 6; ++i) {
    for (std::size_t j = 0; j < 6; ++j) {
      m(i, j) = 10. * std::min(i, j) + std::max(i, j);
    }
  }

  edm4hep::TrackState trackState{};
  EDM4hepUtil::detail::packCovariance(m, trackState.covMatrix);

  for (std::size_t i = 0; i < 6; ++i) {
    for (std::size_t j = 0; j < 6; ++j) {
      BOOST_TEST_INFO_SCOPE("values index (" << i << ", " << j << ")");
      BOOST_CHECK_EQUAL(trackState.getCovMatrix(params.at(i), params.at(j)),
                        static_cast<float>(m(i, j)));
    }
  }
}

// The hit lookup has to be invoked in the order the track states are emitted
// (inside-out), not in the outside-in order they are iterated internally: a
// stateful callback would otherwise associate its hits with the wrong states.
BOOST_AUTO_TEST_CASE(HitLookupInvocationOrder) {
  TrackContainer tracks(std::make_shared<VectorTrackContainer>(),
                        std::make_shared<VectorMultiTrajectory>());
  auto track = tracks.makeTrack();

  BoundVector par;
  par << 1_mm, 5_mm, 0.1, std::numbers::pi / 2. * 0.9, -1 / 1_GeV, 5_ns;
  track.parameters() = par;
  track.covariance() = BoundMatrix::Identity();
  track.setReferenceSurface(
      Surface::makeShared<PerigeeSurface>(Vector3{0, 0, 0}));

  // Appended innermost first, so the internal reverse iteration visits them
  // outermost first.
  const std::array<double, 3> zs{100, 200, 300};
  for (double z : zs) {
    auto ts = track.appendTrackState(TrackStatePropMask::Smoothed);
    ts.typeFlags().setIsMeasurement();
    ts.smoothed() = par;
    ts.smoothedCovariance() = BoundMatrix::Identity();
    ts.setReferenceSurface(
        Surface::makeShared<PerigeeSurface>(Vector3{0, 0, z}));
  }

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext mctx{};
  ConstantBField field{Vector3{0, 0, 2_T}};

  edm4hep::TrackerHit3DCollection hitCollection;
  std::vector<edm4hep::TrackerHit3D> hitByZ;
  for (double z : zs) {
    auto hit = hitCollection.create();
    hit.setPosition({0, 0, z});
    hitByZ.push_back(hit);
  }

  std::vector<double> seen;
  EDM4hepUtil::TrackerHitLookup lookup =
      [&](const AnyConstTrackStateProxy& state)
      -> std::optional<edm4hep::TrackerHit> {
    double z = state.referenceSurface().center(gctx).z();
    seen.push_back(z);
    for (std::size_t i = 0; i < zs.size(); ++i) {
      if (std::abs(z - zs[i]) < 1e-6) {
        return hitByZ[i];
      }
    }
    return std::nullopt;
  };

  edm4hep::TrackCollection outTracks;
  auto to = outTracks.create();
  EDM4hepUtil::writeTrack(gctx, mctx, track, to, field, lookup);

  // Invoked once per measurement state, inside-out
  BOOST_REQUIRE_EQUAL(seen.size(), zs.size());
  for (std::size_t i = 0; i < zs.size(); ++i) {
    BOOST_TEST_INFO_SCOPE("lookup call #" << i);
    CHECK_CLOSE_ABS(seen[i], zs[i], 1e-6);
  }

  // ... and the hits are attached in that same order
  BOOST_REQUIRE_EQUAL(to.trackerHits_size(), zs.size());
  for (std::size_t i = 0; i < zs.size(); ++i) {
    BOOST_TEST_INFO_SCOPE("tracker hit #" << i);
    CHECK_CLOSE_ABS(to.getTrackerHits(i).getPosition().z, zs[i], 1e-6);
  }
}

namespace {
/// Bz grows with |z|, so the field at a track's point of closest approach
/// differs from the field at its perigee reference point.
class ZGradientField : public MagneticFieldProvider {
 public:
  explicit ZGradientField(double b0) : m_b0(b0) {}

  Cache makeCache(const MagneticFieldContext& /*mctx*/) const override {
    return Cache(std::in_place_type<int>, 0);
  }

  Result<Vector3> getField(const Vector3& position,
                           Cache& /*cache*/) const override {
    return Result<Vector3>::success(
        Vector3{0, 0, m_b0 * (1 + 0.001 * std::abs(position.z()))});
  }

 private:
  double m_b0;
};
}  // namespace

// The field is sampled at the point of closest approach, which for parameters
// already on a perigee is d0/z0 away from the reference point. Both sides have
// to sample at the same place, or q/p comes back biased in a non-uniform field.
BOOST_AUTO_TEST_CASE(NonUniformFieldRoundTrip) {
  auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext mctx{};
  ZGradientField field{2_T};

  BoundVector par;
  par << 50_mm, 300_mm, 0.7, 1.1, 1. / 2_GeV, 0;
  BoundMatrix cov = BoundMatrix::Identity() * 1e-6;

  TrackContainer tracks(std::make_shared<VectorTrackContainer>(),
                        std::make_shared<VectorMultiTrajectory>());
  auto track = tracks.makeTrack();
  track.parameters() = par;
  track.covariance() = cov;
  track.setReferenceSurface(
      Surface::makeShared<PerigeeSurface>(Vector3{0, 0, 0}));

  auto ts = track.appendTrackState(TrackStatePropMask::Smoothed);
  ts.typeFlags().setIsMeasurement();
  ts.smoothed() = par;
  ts.smoothedCovariance() = cov;
  ts.setReferenceSurface(
      Surface::makeShared<PerigeeSurface>(Vector3{0, 0, 300_mm}));

  edm4hep::TrackCollection edm4hepTracks;
  auto to = edm4hepTracks.create();
  EDM4hepUtil::writeTrack(gctx, mctx, track, to, field);

  TrackContainer readTracks(std::make_shared<VectorTrackContainer>(),
                            std::make_shared<VectorMultiTrajectory>());
  auto read = readTracks.makeTrack();
  EDM4hepUtil::readTrack(gctx, mctx, to, read, field);

  CHECK_CLOSE_OR_SMALL(track.parameters(), read.parameters(), 1e-5, 1e-8);
}

BOOST_AUTO_TEST_CASE(RoundTripTests) {
  auto trackContainer = std::make_shared<VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  std::mt19937 rng{42};
  std::normal_distribution<double> gauss(0., 1.);
  std::uniform_real_distribution<double> f(-1, 1);
  std::uniform_real_distribution<double> r(0, 1);
  std::uniform_int_distribution<std::uint32_t> nTracks(2, 20);
  std::uniform_int_distribution<std::uint32_t> nTs(1, 20);
  std::uniform_real_distribution<double> phiDist(-std::numbers::pi,
                                                 std::numbers::pi);
  std::uniform_real_distribution<double> etaDist(-4, 4);
  std::uniform_real_distribution<double> ptDist(1_MeV, 10_GeV);
  std::uniform_real_distribution<double> qDist(0., 1.);

  auto genParams = [&]() -> std::pair<BoundVector, BoundMatrix> {
    double d0 = 20_um * gauss(rng);
    double z0 = 20_mm * gauss(rng);
    double phi = phiDist(rng);
    double eta = etaDist(rng);
    double theta = 2 * std::atan(exp(-eta));
    double pt = ptDist(rng);
    double p = pt / std::sin(theta);
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

  std::uint32_t numT = nTracks(rng);
  for (std::uint32_t t = 0; t < numT; t++) {
    auto track = tracks.makeTrack();
    {
      auto [par, cov] = genParams();
      track.parameters() = par;
      track.covariance() = cov;
    }
    track.setReferenceSurface(
        Surface::makeShared<PerigeeSurface>(Vector3{0, 0, 0}));

    std::uint32_t numTs = nTs(rng);
    for (std::uint32_t i = 0; i < numTs; i++) {
      auto ts = track.appendTrackState(TrackStatePropMask::Smoothed);
      double crit = r(rng);
      if (crit < 0.1) {
        ts.typeFlags().setIsHole();
      } else if (crit < 0.2) {
        ts.typeFlags().setIsOutlier();
      } else if (crit < 0.3) {
        ts.typeFlags().setIsSharedHit();
      } else if (crit < 0.4) {
        ts.typeFlags().setIsMaterial();
      } else {
        ts.typeFlags().setIsMeasurement();
      }

      auto [par, cov] = genParams();
      ts.smoothed() = par;
      ts.smoothedCovariance() = cov;
      Vector3 pos;
      pos << 1000 * f(rng), 1000 * f(rng), 3000 * f(rng);
      ts.setReferenceSurface(Surface::makeShared<PerigeeSurface>(pos));
    }

    calculateTrackQuantities(track);
  }

  edm4hep::TrackCollection edm4hepTracks;

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  double Bz = 3_T;

  auto logger = getDefaultLogger("EDM4hep", Logging::INFO);

  for (const auto& track : tracks) {
    auto to = edm4hepTracks.create();
    EDM4hepUtil::writeTrack(gctx, track, to, Bz, {}, *logger);
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

  TrackContainer readTracks(std::make_shared<VectorTrackContainer>(),
                            std::make_shared<VectorMultiTrajectory>());

  for (const auto edm4hepTrack : edm4hepTracksConst) {
    auto track = readTracks.makeTrack();
    EDM4hepUtil::readTrack(gctx, edm4hepTrack, track, Bz, *logger);
  }

  BOOST_CHECK_EQUAL(tracks.size(), readTracks.size());
  std::size_t t = 0;

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

    // Track summary quantities have to survive the round trip as well
    CHECK_CLOSE_OR_SMALL(orig.chi2(), read.chi2(), 1e-5, 1e-8);
    BOOST_CHECK_EQUAL(orig.nDoF(), read.nDoF());
    BOOST_CHECK_EQUAL(orig.nHoles(), read.nHoles());

    auto origTsIt = orig.trackStatesReversed().begin();
    auto readTsIt = read.trackStatesReversed().begin();

    std::size_t tsi = 0;
    while (origTsIt != orig.trackStatesReversed().end() &&
           readTsIt != read.trackStatesReversed().end()) {
      BOOST_TEST_INFO_SCOPE("TS: #" << tsi);
      auto nextMeas = std::find_if(
          origTsIt, orig.trackStatesReversed().end(),
          [](const auto& ts) { return ts.typeFlags().isMeasurement(); });
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

}  // namespace ActsTests
