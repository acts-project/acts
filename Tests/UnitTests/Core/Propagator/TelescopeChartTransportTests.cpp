// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/DefaultActor.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <cmath>

namespace ActsTests {
namespace {

using namespace Acts;

const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext mctx;
const auto telescope = BoundParametrization::telescope();
const auto canonical = BoundParametrization::globalSpherical();
const auto pion = ParticleHypothesis::pion();
using Prop = Propagator<StraightLineStepper, VoidNavigator>;

auto plane(double z, const BoundParametrization& chart) {
  Transform3 transform = Transform3::Identity();
  transform.translation().z() = z;
  return Surface::makeShared<PlaneSurface>(transform, nullptr, chart);
}

BoundMatrix seedCovariance() {
  BoundMatrix covariance = BoundMatrix::Identity();
  covariance(eBoundLoc0, eBoundDir1) = covariance(eBoundDir1, eBoundLoc0) =
      0.25;
  covariance(eBoundLoc1, eBoundDir0) = covariance(eBoundDir0, eBoundLoc1) =
      -0.125;
  return covariance;
}

// Independent map between two parallel xy planes, including flight time for a
// singly charged pion. Both frames are fixed under translation along z.
BoundVector intersect(const BoundVector& start, double dz, bool sourceRotated,
                      bool targetRotated) {
  const double phi = start[eBoundDir0];
  const double theta = start[eBoundDir1];
  Vector3 d(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
            std::cos(theta));
  if (sourceRotated) {
    d = Vector3(d.z(), d.x(), d.y());
  }
  const double length = dz / d.z();
  BoundVector end = start;
  end.head<2>() += length * d.head<2>();
  end[eBoundTime] += length * std::hypot(1., pion.mass() * start[eBoundQOverP]);
  if (targetRotated) {
    d = Vector3(d.y(), d.z(), d.x());
  }
  end[eBoundDir0] = std::atan2(d.y(), d.x());
  end[eBoundDir1] = std::atan2(std::hypot(d.x(), d.y()), d.z());
  return end;
}

// Compare the actual stepper and covariance transport with finite differences
// of the analytic map (including intersection, time, and both chart bases).
auto checkTransport(const BoundTrackParameters& start,
                    const std::shared_ptr<PlaneSurface>& target) {
  const auto& src = start.referenceSurface();
  const double dz = target->center(gctx).z() - src.center(gctx).z();
  const auto map = [&](const BoundVector& v) {
    return intersect(v, dz, src.boundParametrization(gctx) == telescope,
                     target->boundParametrization(gctx) == telescope);
  };
  const auto& params = start.parameters();
  BoundMatrix numerical;
  constexpr double h = 1e-5;
  for (Eigen::Index i = 0; i < params.size(); ++i) {
    const BoundVector delta = h * BoundVector::Unit(i);
    numerical.col(i) = (map(params + delta) - map(params - delta)) / (2. * h);
  }

  Prop prop{StraightLineStepper{}, VoidNavigator{}};
  Prop::Options<> options(gctx, mctx);
  options.direction =
      Direction::fromScalarZeroAsPositive(dz / start.direction(gctx).z());
  auto result = prop.propagate(start, *target, options);
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE(result->transportJacobian && result->endParameters);
  CHECK_CLOSE_ABS(*result->transportJacobian, numerical, 2e-7);
  const auto& end = *result->endParameters;
  CHECK_CLOSE_ABS(end.parameters(), map(params), 1e-9);
  BOOST_REQUIRE(end.covariance() && end.covariance()->allFinite());
  CHECK_CLOSE_COVARIANCE(
      *end.covariance(),
      numerical * *start.covariance() * numerical.transpose(), 1e-8);
  return end;
}

}  // namespace

BOOST_AUTO_TEST_CASE(TelescopeTransportThroughCurvilinearAnchor) {
  const auto source = plane(0., telescope);
  Prop prop{StraightLineStepper{}, VoidNavigator{}};
  for (double sign : {-1., 1.}) {
    for (double tilt : {0., 1e-5}) {
      const Vector3 d = Vector3(tilt, -2. * tilt, sign).normalized();
      auto seed = BoundTrackParameters::create(gctx, source, Vector4::Zero(), d,
                                               0.5, seedCovariance(), pion);
      BOOST_REQUIRE(seed.ok());
      for (auto direction : {Direction::Forward(), Direction::Backward()}) {
        const auto target = plane(30. * sign * direction, telescope);
        const auto direct = checkTransport(*seed, target);
        Prop::Options<> options(gctx, mctx);
        options.direction = direction;
        options.pathLimit = 10.;
        // With no target the propagator transports covariance to a curvilinear
        // anchor and returns bound parameters on a synthesized plane.
        const auto mid = prop.propagate(*seed, options);
        BOOST_REQUIRE(mid.ok());
        BOOST_REQUIRE(mid->endParameters);
        const auto& anchor = *mid->endParameters;
        BOOST_REQUIRE(anchor.referenceSurface().boundParametrization(gctx) ==
                      telescope);
        CHECK_CLOSE_ABS(anchor.direction(gctx), d, 1e-12);
        options.pathLimit = 50.;
        const auto via = prop.propagate(anchor, *target, options);
        BOOST_REQUIRE(via.ok());
        BOOST_REQUIRE(via->endParameters && via->endParameters->covariance());
        BOOST_REQUIRE(via->endParameters->covariance()->allFinite());
        CHECK_CLOSE_ABS(via->endParameters->parameters(), direct.parameters(),
                        1e-9);
        CHECK_CLOSE_COVARIANCE(*via->endParameters->covariance(),
                               *direct.covariance(), 1e-9);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ScatterInCurvilinearTelescopeChart) {
  using MaterialPropagator = Propagator<StraightLineStepper, DirectNavigator>;
  MaterialPropagator propagator{StraightLineStepper{}, DirectNavigator{}};
  const auto start = plane(0., telescope);
  // The passive slab's chart differs from that of the active covariance.
  auto slab = plane(10., canonical);
  slab->assignSurfaceMaterial(std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(makeSilicon(), 0.3)));
  const auto target = plane(30., telescope);
  const auto seed = BoundTrackParameters::create(gctx, start, Vector4::Zero(),
                                                 Vector3::UnitZ(), 0.5,
                                                 seedCovariance(), pion);
  BOOST_REQUIRE(seed.ok());
  using Actors = ActorList<DefaultActor, MaterialInteractor>;
  MaterialPropagator::Options<Actors> options(gctx, mctx);
  options.navigation.surfaceSequence = {slab.get(), target.get()};
  auto& material = options.actor.get<MaterialInteractor>();
  material.energyLoss = false;
  material.multipleScattering = false;
  material.recordInteractions = true;
  const auto vacuum = propagator.propagate(*seed, *target, options);
  BOOST_REQUIRE(vacuum.ok() && vacuum->endParameters &&
                vacuum->endParameters->covariance());
  material.multipleScattering = true;
  auto state = propagator.makeState(options);
  const auto result = propagator.propagate(state, *seed, target.get());
  BOOST_REQUIRE(result.ok() && result->endParameters &&
                result->endParameters->covariance());
  const auto& interactions =
      state.actor.get<MaterialInteractor::State>().materialInteractions;
  BOOST_REQUIRE_EQUAL(interactions.size(), 1u);
  const auto& hit = interactions.front();
  BOOST_CHECK(std::isinf(hit.sigmaPhi2));  // legacy global-pole diagnostic
  BOOST_REQUIRE_GT(hit.sigmaTheta2, 0.);
  const BoundMatrix added = *result->endParameters->covariance() -
                            *vacuum->endParameters->covariance();
  BOOST_CHECK(added.allFinite());
  CHECK_CLOSE_ABS((added.block<2, 2>(eBoundDir0, eBoundDir0)),
                  SquareMatrix2::Identity() * hit.sigmaTheta2, 1e-12);
}

BOOST_AUTO_TEST_CASE(MixedTelescopeTransport) {
  // In this direction both charts are regular and have distinct local bases.
  auto seed = BoundTrackParameters::create(
      gctx, plane(0., canonical), Vector4::Zero(),
      Vector3(0.3, 0.2, 1.).normalized(), 0.5, seedCovariance(), pion);
  BOOST_REQUIRE(seed.ok());
  const auto end = plane(30., canonical);
  const auto direct = checkTransport(*seed, end);
  const auto middle = checkTransport(*seed, plane(10., telescope));
  const auto mixed = checkTransport(middle, end);
  CHECK_CLOSE_ABS(mixed.parameters(), direct.parameters(), 1e-10);
  CHECK_CLOSE_COVARIANCE(*mixed.covariance(), *direct.covariance(), 1e-9);
}

}  // namespace ActsTests
