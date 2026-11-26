// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldError.hpp"
#include "Acts/MagneticField/MagneticFieldQuery.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/TwoTrackVertexFinder.hpp"
#include "Acts/Vertexing/TwoTrackVertexSeed.hpp"
#include "Acts/Vertexing/VertexingError.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numbers>
#include <random>
#include <stdexcept>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {
namespace {
const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext mctx;

BoundTrackParameters trackAt(const Vector3& position, const Vector3& direction,
                             double qOverP) {
  auto surface = Surface::makeShared<PerigeeSurface>(position);
  Vector4 pos;
  pos.head<3>() = position;
  pos[eTime] = 0.;
  return BoundTrackParameters::create(
             gctx, surface, pos, direction.normalized(), qOverP,
             BoundMatrix::Identity(), ParticleHypothesis::pion())
      .value();
}

class TestFitter final : public IVertexFitter {
 public:
  mutable unsigned int calls = 0;
  bool fail = false;
  Result<Vertex> fit(const VertexFitInput& input, const GeometryContext& geo,
                     MagneticFieldQuery& field) const override {
    ++calls;
    BOOST_CHECK_EQUAL(&geo, &gctx);
    BOOST_CHECK_EQUAL(&field.context(), &mctx);
    BOOST_CHECK_EQUAL(input.tracks.size(), 2u);
    BOOST_CHECK(!input.constraint.has_value());
    if (fail) {
      return VertexingError::NotConverged;
    }
    Vertex vertex(input.seedPosition);
    vertex.setFullCovariance(SquareMatrix4::Identity());
    vertex.setFitQuality({1., 1.});
    return vertex;
  }
};

class FailingField final : public MagneticFieldProvider {
 public:
  Cache makeCache(const MagneticFieldContext& /*mctx*/) const override {
    return Cache(std::in_place_type<int>, 0);
  }
  Result<Vector3> getField(const Vector3& /*position*/,
                           Cache& /*cache*/) const override {
    return Result<Vector3>::failure(MagneticFieldError::OutOfBounds);
  }
};
}  // namespace

BOOST_AUTO_TEST_SUITE(VertexingSuite)

BOOST_AUTO_TEST_CASE(two_track_straight_seed) {
  const Vector3 truth(20., -5., 12.);
  const auto a = trackAt(truth - Vector3(2., 0., 0.), Vector3::UnitX(), 1.);
  const auto b = trackAt(truth - Vector3(0., 3., 0.), Vector3::UnitY(), -1.);
  auto seed = estimateTwoTrackVertexSeed(a, b, gctx, 0.);
  BOOST_REQUIRE(seed.ok());
  CHECK_CLOSE_ABS(seed->position, truth, 1e-12);
  BOOST_CHECK_SMALL(seed->separation, 1e-12);
  const auto skew =
      trackAt(truth + Vector3(0., -3., 2.), Vector3::UnitY(), -1.);
  seed = estimateTwoTrackVertexSeed(a, skew, gctx, 0.);
  BOOST_REQUIRE(seed.ok());
  CHECK_CLOSE_ABS(seed->position, truth + Vector3(0., 0., 1.), 1e-12);
  CHECK_CLOSE_ABS(seed->separation, 2., 1e-12);
  const auto parallel = trackAt(Vector3::Zero(), Vector3(1., 1e-8, 0.), -1.);
  auto failure = estimateTwoTrackVertexSeed(a, parallel, gctx, 0.);
  BOOST_REQUIRE(!failure.ok());
  BOOST_CHECK(failure.error() == VertexingError::SeedingError);
  BOOST_CHECK(!estimateTwoTrackVertexSeed(
                   a, b, gctx, std::numeric_limits<double>::quiet_NaN())
                   .ok());
}

BOOST_AUTO_TEST_CASE(two_track_circle_geometry) {
  // Unit axial field and q/p=-1/r give counterclockwise circles.
  const auto first = trackAt(Vector3(5., 0., 0.), Vector3::UnitY(), -0.2);
  const auto second = trackAt(Vector3(13., 0., 0.), Vector3::UnitY(), -0.2);
  auto seed = estimateTwoTrackVertexSeed(first, second, gctx, 1.);
  BOOST_REQUIRE(seed.ok());
  CHECK_CLOSE_ABS(seed->position.x(), 4., 1e-12);
  CHECK_CLOSE_ABS(std::abs(seed->position.y()), 3., 1e-12);
  BOOST_CHECK_SMALL(seed->separation, 1e-12);
  const auto reversed = estimateTwoTrackVertexSeed(second, first, gctx, 1.);
  BOOST_REQUIRE(reversed.ok());
  CHECK_CLOSE_ABS(reversed->separation, seed->separation, 1e-12);
  // Reverse both charges and the field: the trajectories are unchanged.
  const auto positive1 = trackAt(Vector3(5., 0., 0.), Vector3::UnitY(), 0.2);
  const auto positive2 = trackAt(Vector3(13., 0., 0.), Vector3::UnitY(), 0.2);
  const auto flipped =
      estimateTwoTrackVertexSeed(positive1, positive2, gctx, -1.);
  BOOST_REQUIRE(flipped.ok());
  CHECK_CLOSE_ABS(flipped->position, seed->position, 1e-12);
  const auto disjoint = trackAt(Vector3(17., 0., 0.), Vector3::UnitY(), -0.2);
  seed = estimateTwoTrackVertexSeed(first, disjoint, gctx, 1.);
  BOOST_REQUIRE(seed.ok());
  CHECK_CLOSE_ABS(seed->position, Vector3(6., 0., 0.), 1e-12);
  CHECK_CLOSE_ABS(seed->separation, 2., 1e-12);
  const auto nested = trackAt(Vector3(3., 0., 0.), Vector3::UnitY(), -0.5);
  seed = estimateTwoTrackVertexSeed(first, nested, gctx, 1.);
  BOOST_REQUIRE(seed.ok());
  CHECK_CLOSE_ABS(seed->position, Vector3(4., 0., 0.), 1e-12);
  CHECK_CLOSE_ABS(seed->separation, 2., 1e-12);
  const auto tangent = trackAt(Vector3(15., 0., 0.), Vector3::UnitY(), -0.2);
  seed = estimateTwoTrackVertexSeed(first, tangent, gctx, 1.);
  BOOST_REQUIRE(seed.ok());
  CHECK_CLOSE_ABS(seed->position, Vector3(5., 0., 0.), 1e-12);
  BOOST_CHECK_SMALL(seed->separation, 1e-12);
  auto failure = estimateTwoTrackVertexSeed(first, first, gctx, 1.);
  BOOST_REQUIRE(!failure.ok());
  BOOST_CHECK(failure.error() == VertexingError::SeedingError);
}

BOOST_AUTO_TEST_CASE(two_track_finder_contract) {
  auto fitter = std::make_shared<TestFitter>();
  auto field = std::make_shared<ConstantBField>(Vector3::Zero());
  TwoTrackVertexFinder::Config cfg;
  cfg.vertexFitter = fitter;
  cfg.field = field;
  cfg.extractParameters.connect<&InputTrack::extractParameters>();
  auto invalid = cfg;
  invalid.vertexFitter.reset();
  BOOST_CHECK_THROW(TwoTrackVertexFinder{invalid}, std::invalid_argument);
  TwoTrackVertexFinder finder(cfg);
  const IVertexFinder& iface = finder;
  auto state = iface.makeState(mctx);
  const VertexingOptions options(gctx, mctx);
  const Vector3 truth(20., -5., 12.);
  const auto first = trackAt(truth, Vector3::UnitX(), 1.);
  const auto second = trackAt(truth, Vector3::UnitY(), -1.);
  const std::vector<InputTrack> tracks{InputTrack(&first), InputTrack(&second)};
  BOOST_CHECK(iface.find({}, options, state)->empty());
  auto result = iface.find(tracks, options, state);
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 1u);
  CHECK_CLOSE_ABS(result->front().position(), truth, 1e-12);
  BOOST_CHECK_EQUAL(fitter->calls, 1u);
  auto duplicate = tracks;
  duplicate.push_back(tracks[0]);
  result = iface.find(duplicate, options, state);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result->size(), 1u);
  iface.setTracksToRemove(state, {tracks[0]});
  BOOST_CHECK(iface.find(tracks, options, state)->empty());
  state = iface.makeState(mctx);
  fitter->fail = true;
  result = iface.find(tracks, options, state);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == VertexingError::NotConverged);
  fitter->fail = false;
  cfg.field = std::make_shared<FailingField>();
  TwoTrackVertexFinder badFieldFinder(cfg);
  state = badFieldFinder.makeState(mctx);
  result = badFieldFinder.find(tracks, options, state);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == MagneticFieldError::OutOfBounds);
  cfg.field = field;
  cfg.maxChi2 = 0.;
  TwoTrackVertexFinder tightFinder(cfg);
  state = tightFinder.makeState(mctx);
  BOOST_CHECK(tightFinder.find(tracks, options, state)->empty());
  Vertex constraint(truth);
  constraint.setFullCovariance(SquareMatrix4::Identity());
  const VertexingOptions constrained(gctx, mctx, constraint);
  state = iface.makeState(mctx);
  result = iface.find(tracks, constrained, state);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == VertexingError::InvalidInput);
  cfg.maxChi2 = 25.;
  cfg.field = std::make_shared<ConstantBField>(Vector3(1_T, 0., 0.));
  TwoTrackVertexFinder tiltedFieldFinder(cfg);
  state = tiltedFieldFinder.makeState(mctx);
  result = tiltedFieldFinder.find(tracks, options, state);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == VertexingError::InvalidInput);
  const auto sameCharge = trackAt(truth, Vector3::UnitY(), 1.);
  state = iface.makeState(mctx);
  BOOST_CHECK(iface.find({tracks[0], InputTrack(&sameCharge)}, options, state)
                  ->empty());
  const auto distant =
      trackAt(truth + Vector3(0., 0., 10.), Vector3::UnitY(), -1.);
  BOOST_CHECK(
      iface.find({tracks[0], InputTrack(&distant)}, options, state)->empty());
  const auto background =
      trackAt(truth + Vector3(0., 0., 20.), Vector3::UnitX(), 1.);
  result = iface.find(
      {tracks[0], tracks[1], InputTrack(&distant), InputTrack(&background)},
      options, state);
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 1u);
  CHECK_CLOSE_ABS(result->front().position(), truth, 1e-12);
}

BOOST_AUTO_TEST_CASE(displaced_two_body_pion_gun) {
  // Correlated pion daughters only: no parent transport or decay simulation.
  constexpr double parentMass = 0.497611_GeV;
  constexpr double pionMass = 0.13957039_GeV;
  const double restEnergy = parentMass / 2.;
  const double restMomentum =
      std::sqrt(restEnergy * restEnergy - pionMass * pionMass);
  auto field = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));
  auto propagator =
      std::make_shared<Propagator<EigenStepper<>>>(EigenStepper<>(field));
  HelicalTrackLinearizer::Config linearizerCfg;
  linearizerCfg.propagator = propagator;
  HelicalTrackLinearizer linearizer(linearizerCfg);
  FullBilloirVertexFitter::Config fitCfg;
  fitCfg.extractParameters.connect<&InputTrack::extractParameters>();
  fitCfg.trackLinearizer.connect<&HelicalTrackLinearizer::linearizeTrack>(
      &linearizer);
  TwoTrackVertexFinder::Config cfg;
  cfg.vertexFitter = std::make_shared<FullBilloirVertexFitter>(fitCfg);
  cfg.field = field;
  cfg.extractParameters.connect<&InputTrack::extractParameters>();
  TwoTrackVertexFinder finder(cfg);
  const IVertexFinder& iface = finder;
  const VertexingOptions options(gctx, mctx);
  auto origin = Surface::makeShared<PerigeeSurface>(Vector3::Zero());
  std::mt19937 rng(314159);
  std::normal_distribution<double> gaussian;
  std::uniform_real_distribution<double> uniform(-1., 1.);
  unsigned int accepted = 0;
  unsigned int total = 0;
  double sumMass = 0.;
  Vector3 sumPull = Vector3::Zero();
  Vector3 sumPullSquared = Vector3::Zero();
  double sumResidualSquared = 0.;
  for (double radius : {1_mm, 10_mm, 50_mm, 100_mm}) {
    for (double momentum : {1_GeV, 2_GeV, 5_GeV}) {
      unsigned int binAccepted = 0;
      for (unsigned int event = 0; event < 20; ++event) {
        const double azimuth = std::numbers::pi * uniform(rng);
        const Vector3 truth(radius * std::cos(azimuth),
                            radius * std::sin(azimuth), 20_mm * uniform(rng));
        const Vector3 axis(std::cos(azimuth), std::sin(azimuth), 0.);
        const double energy = std::hypot(momentum, parentMass);
        const double gamma = energy / parentMass;
        const double beta = momentum / energy;
        const double cosTheta = 0.8 * uniform(rng);
        const double phi = std::numbers::pi * uniform(rng);
        const Vector3 rest =
            restMomentum *
            Vector3(std::sqrt(1. - cosTheta * cosTheta) * std::cos(phi),
                    std::sqrt(1. - cosTheta * cosTheta) * std::sin(phi),
                    cosTheta);
        std::vector<BoundTrackParameters> tracks;
        std::vector<BoundTrackParameters> exactTracks;
        for (double charge : {1., -1.}) {
          const Vector3 daughter = charge * rest;
          const Vector3 boosted =
              daughter +
              ((gamma - 1.) * daughter.dot(axis) + gamma * beta * restEnergy) *
                  axis;
          auto atVertex = trackAt(truth, boosted, charge / boosted.norm());
          PropagatorPlainOptions propagation(gctx, mctx);
          propagation.direction = Direction::Backward();
          auto atOrigin =
              propagator->propagateToSurface(atVertex, *origin, propagation);
          BOOST_REQUIRE(atOrigin.ok());
          exactTracks.push_back(*atOrigin);
          BoundVector measured = atOrigin->parameters();
          BoundVector sigma;
          sigma << 0.02_mm, 0.05_mm, 0.0002, 0.0002,
              std::abs(measured[eBoundQOverP]) * 0.001, 0.1_ns;
          BoundMatrix covariance = sigma.cwiseAbs2().asDiagonal();
          for (unsigned int i = 0; i < eBoundSize; ++i) {
            measured[i] += sigma[i] * gaussian(rng);
          }
          tracks.emplace_back(origin, measured, covariance,
                              ParticleHypothesis::pion());
        }
        const auto exactSeed = estimateTwoTrackVertexSeed(
            exactTracks[0], exactTracks[1], gctx, 2_T);
        BOOST_REQUIRE(exactSeed.ok());
        CHECK_CLOSE_ABS(exactSeed->position, truth, 1e-3_mm);
        BOOST_CHECK_LT(exactSeed->separation, 1e-3_mm);
        const std::vector<InputTrack> input{InputTrack(&tracks[0]),
                                            InputTrack(&tracks[1])};
        auto state = iface.makeState(mctx);
        auto result = iface.find(input, options, state);
        BOOST_REQUIRE(result.ok());
        ++total;
        if (result->empty()) {
          continue;
        }
        BOOST_REQUIRE_EQUAL(result->size(), 1u);
        const auto& vertex = result->front();
        BOOST_REQUIRE_EQUAL(vertex.tracks().size(), 2u);
        ++accepted;
        ++binAccepted;
        const Vector3 residual = vertex.position() - truth;
        sumResidualSquared += residual.squaredNorm();
        BOOST_REQUIRE((vertex.covariance().diagonal().array() > 0.).all());
        const Vector3 pull =
            residual.cwiseQuotient(vertex.covariance().diagonal().cwiseSqrt());
        sumPull += pull;
        sumPullSquared += pull.cwiseAbs2();
        Vector3 totalMomentum = Vector3::Zero();
        double totalEnergy = 0.;
        for (const auto& track : vertex.tracks()) {
          const Vector3 fittedMomentum = track.fittedParams.momentum();
          totalMomentum += fittedMomentum;
          totalEnergy += std::hypot(fittedMomentum.norm(), pionMass);
        }
        const double mass = std::sqrt(std::max(
            0., totalEnergy * totalEnergy - totalMomentum.squaredNorm()));
        sumMass += mass;
        CHECK_CLOSE_ABS(mass, parentMass, 0.02_GeV);
      }
      BOOST_CHECK_GE(binAccepted, 16u);
      BOOST_TEST_MESSAGE("radius=" << radius / 1_mm
                                   << " mm, p=" << momentum / 1_GeV << " GeV: "
                                   << binAccepted << "/20 accepted");
    }
  }
  BOOST_REQUIRE_GT(accepted, 0u);
  const double rms = std::sqrt(sumResidualSquared / accepted);
  const Vector3 meanPull = sumPull / accepted;
  const Vector3 rmsPull = (sumPullSquared / accepted).cwiseSqrt();
  BOOST_TEST_MESSAGE("Displaced pion gun: "
                     << accepted << "/" << total
                     << ", vertex RMS=" << rms / 1_mm
                     << " mm, mean mass=" << sumMass / accepted / 1_GeV
                     << " GeV, pull means=" << meanPull.transpose()
                     << ", pull RMS=" << rmsPull.transpose());
  BOOST_CHECK_LT(rms, 2_mm);
  CHECK_CLOSE_ABS(sumMass / accepted, parentMass, 0.003_GeV);
  for (unsigned int i = 0; i < 3; ++i) {
    BOOST_CHECK_LT(std::abs(meanPull[i]), 0.4);
    BOOST_CHECK_GT(rmsPull[i], 0.5);
    BOOST_CHECK_LT(rmsPull[i], 1.5);
  }
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
