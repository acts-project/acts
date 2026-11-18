// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/PointSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/detail/CurvatureCorrection.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <limits>

#include "SurfaceStub.hpp"

using namespace Acts;
namespace ActsTests {
const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

BOOST_AUTO_TEST_CASE(CurvatureCorrectionModes) {
  const Vector3 direction = Vector3::UnitZ();
  const Vector3 normal = Vector3(1., 0., 1.).normalized();
  const double s = 100.;
  const Vector3 seed = s * direction;
  for (double k : {0., 0.0005, -0.0005}) {
    const Vector3 curvature(k, 0., 0.);
    for (double scale : {1., 100.}) {
      const auto first = detail::correctCurvatureIntersection<false>(
          s, seed, direction, curvature, scale * normal.dot(direction),
          scale * normal.dot(curvature), curvature.norm(), scale);
      const auto quadratic = detail::correctCurvatureIntersection<true>(
          s, seed, direction, curvature, scale * normal.dot(direction),
          scale * normal.dot(curvature), curvature.norm(), scale);
      BOOST_REQUIRE(first);
      BOOST_REQUIRE(quadratic);
      CHECK_CLOSE_ABS(normal.dot(first->position - seed), 0., 1e-12);
      CHECK_CLOSE_ABS(normal.dot(quadratic->position - seed), 0., 1e-12);
      const double path = quadratic->pathLength;
      CHECK_CLOSE_ABS(quadratic->position,
                      (path * direction + 0.5 * path * path * curvature).eval(),
                      1e-12);
      const auto selected = detail::correctCurvatureIntersection(
          s, seed, direction, curvature, scale * normal.dot(direction),
          scale * normal.dot(curvature), curvature.norm(), scale);
      BOOST_REQUIRE(selected);
      const auto& expected =
          detail::s_useQuadraticCurvatureCorrection ? quadratic : first;
      CHECK_CLOSE_ABS(selected->position, expected->position, 1e-12);
      CHECK_CLOSE_ABS(selected->pathLength, expected->pathLength, 1e-12);
    }
  }
  // A local parabola that turns before reaching the plane has no real root.
  BOOST_CHECK(!detail::correctCurvatureIntersection<true>(
      s, seed, direction, Vector3(-0.001, 0., 0.), 0.1, -0.001, 0.001));
}

// Compare against an independent constant-field circle, not the truncated
// trajectory. The quadratic solve need not improve the missing helix terms.
BOOST_AUTO_TEST_CASE(CurvatureQuadraticComparison) {
  const Vector3 start = Vector3::Zero();
  const Vector3 direction = Vector3::UnitZ();
  for (double tilt : {0., 0.5, 1.}) {
    Transform3 transform = Transform3::Identity();
    transform.translation() = Vector3(0., 0., 200.);
    transform.linear() = AngleAxis3(tilt, Vector3::UnitY()).toRotationMatrix();
    const auto plane = Surface::makeShared<PlaneSurface>(
        transform, std::make_shared<RectangleBounds>(1000., 1000.));
    for (double k : {0., 0.00025, 0.0005, -0.0005}) {
      const Vector3 curvature(k, 0., 0.);
      const auto seed = plane->intersect(gctx, start, direction).at(0);
      const Vector3 normal = transform.linear().col(2);
      const auto first =
          detail::PlanarHelper::correctPlanarCurvatureIntersection<false>(
              seed, direction, curvature, normal);
      const auto quadratic =
          detail::PlanarHelper::correctPlanarCurvatureIntersection<true>(
              seed, direction, curvature, normal);
      const auto selected = plane
                                ->intersect(gctx, start, direction, curvature,
                                            BoundaryTolerance::None())
                                .at(0);
      const auto& expected =
          detail::s_useQuadraticCurvatureCorrection ? quadratic : first;
      CHECK_CLOSE_ABS(selected.position(), expected.position(), 1e-10);
      CHECK_CLOSE_ABS(selected.pathLength(), expected.pathLength(), 1e-10);
      BOOST_REQUIRE(first.isValid());
      BOOST_REQUIRE(quadratic.isValid());
      const double path = quadratic.pathLength();
      CHECK_CLOSE_ABS(quadratic.position(),
                      (path * direction + 0.5 * path * path * curvature).eval(),
                      1e-10);
      BOOST_CHECK(plane->isOnSurface(gctx, quadratic.position(), direction));
      const auto reverseSeed = plane->intersect(gctx, start, -direction).at(0);
      const auto reverse =
          detail::PlanarHelper::correctPlanarCurvatureIntersection<true>(
              reverseSeed, -direction, curvature, normal);
      CHECK_CLOSE_ABS(reverse.position(), quadratic.position(), 1e-10);
      CHECK_CLOSE_ABS(reverse.pathLength(), -path, 1e-10);
      const auto exact = [k](double s) -> Vector3 {
        if (k == 0.) {
          return Vector3(0., 0., s);
        }
        return Vector3(2. * std::pow(std::sin(k * s / 2.), 2) / k, 0.,
                       std::sin(k * s) / k);
      };
      double lo = 0.;
      double hi = 400.;
      for (int i = 0; i < 60; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (normal.dot(exact(mid) - transform.translation()) < 0.) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      const Vector3 reference = exact(0.5 * (lo + hi));
      std::cout << "tilt=" << tilt << " k=" << k
                << " first_error_mm=" << (first.position() - reference).norm()
                << " quadratic_error_mm="
                << (quadratic.position() - reference).norm() << '\n';
      if (k == 0.) {
        CHECK_CLOSE_ABS(quadratic.position(), first.position(), 1e-10);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(CurvaturePlaneBoundsAndDistance) {
  const auto bounds = std::make_shared<RectangleBounds>(0.2, 0.2);
  const auto sensor = Surface::makeShared<PlaneSurface>(
      Transform3(Translation3(1., 0., 20.)), bounds);
  const Vector3 start = Vector3::Zero();
  const Vector3 direction = Vector3::UnitZ();
  const Vector3 curvature(0.005, 0., 0.);
  BOOST_CHECK(
      !sensor->intersect(gctx, start, direction, BoundaryTolerance::None())
           .at(0)
           .isValid());
  const auto crossing = sensor->intersect(gctx, start, direction, curvature,
                                          BoundaryTolerance::None());
  BOOST_REQUIRE(crossing.at(0).isValid());
  BOOST_CHECK_EQUAL(crossing.size(), 1u);
  CHECK_CLOSE_ABS(crossing.at(0).position(), Vector3(1., 0., 20.), 1e-12);
  CHECK_CLOSE_ABS(crossing.at(0).pathLength(), 20., 1e-12);
  // Analytic circle with initial tangent +z and curvature towards +x.
  const double exactX = (1. - std::sqrt(1. - std::pow(20. * 0.005, 2))) / 0.005;
  BOOST_CHECK_LT(std::abs(crossing.at(0).position().x() - exactX), 0.003);
  const auto reversed = sensor->intersect(gctx, start, -direction, curvature,
                                          BoundaryTolerance::None());
  CHECK_CLOSE_ABS(reversed.at(0).position(), crossing.at(0).position(), 1e-12);
  CHECK_CLOSE_ABS(reversed.at(0).pathLength(), -20., 1e-12);
  BOOST_CHECK(!sensor
                   ->intersect(gctx, start, direction, -curvature,
                               BoundaryTolerance::None())
                   .at(0)
                   .isValid());

  Transform3 tilted = Transform3::Identity();
  tilted.translation() = Vector3(0., 0., 20.);
  tilted.linear() =
      AngleAxis3(std::atan(0.2), Vector3::UnitY()).toRotationMatrix();
  const auto plane = Surface::makeShared<PlaneSurface>(
      tilted, std::make_shared<RectangleBounds>(10., 10.));
  const auto hit = plane
                       ->intersect(gctx, start, direction, curvature,
                                   BoundaryTolerance::None())
                       .at(0);
  BOOST_REQUIRE(hit.isValid());
  CHECK_CLOSE_ABS(hit.pathLength(), 19.8, 1e-12);
  CHECK_CLOSE_ABS(hit.position(), Vector3(1., 0., 19.8), 1e-12);
  BOOST_CHECK(plane->isOnSurface(gctx, hit.position(), direction));
}

BOOST_AUTO_TEST_CASE(CurvatureCylinderAndFallbacks) {
  const auto cylinder =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20., 10.);
  const Vector3 direction = Vector3::UnitX();
  const auto hits =
      cylinder->intersect(gctx, Vector3::Zero(), direction,
                          Vector3(0., 0.005, 0.), BoundaryTolerance::None());
  BOOST_CHECK_EQUAL(hits.size(), 2u);
  for (const auto& hit : hits) {
    BOOST_REQUIRE(hit.isValid());
    CHECK_CLOSE_ABS(hit.position().head<2>().norm(), 20., 1e-12);
    BOOST_CHECK_GT(hit.position().y(), 0.);
  }
  CHECK_CLOSE_ABS(hits.at(0).pathLength(), -20., 1e-12);
  CHECK_CLOSE_ABS(hits.at(1).pathLength(), 20., 1e-12);
  Transform3 transform = Transform3::Identity();
  transform.translation() = Vector3(3., -4., 7.);
  transform.linear() =
      AngleAxis3(0.7, Vector3(1., 2., 3.).normalized()).toRotationMatrix();
  const auto moved = Surface::makeShared<CylinderSurface>(transform, 20., 10.);
  const auto movedHits = moved->intersect(
      gctx, transform * Vector3::Zero(), transform.linear() * direction,
      transform.linear() * Vector3(0., 0.005, 0.), BoundaryTolerance::None());
  CHECK_CLOSE_ABS(movedHits.at(1).position(), transform * hits.at(1).position(),
                  1e-12);

  const auto plane = Surface::makeShared<PlaneSurface>(
      Transform3(Translation3(0., 0., 20.)),
      std::make_shared<RectangleBounds>(10., 10.));
  for (const Vector3& curvature :
       {Vector3::Zero().eval(), Vector3(1., 0., 0.)}) {
    const auto straight =
        plane->intersect(gctx, Vector3::Zero(), Vector3::UnitZ());
    const auto fallback =
        plane->intersect(gctx, Vector3::Zero(), Vector3::UnitZ(), curvature,
                         BoundaryTolerance::None());
    CHECK_CLOSE_ABS(fallback.at(0).position(), straight.at(0).position(),
                    1e-12);
    BOOST_CHECK_EQUAL(fallback.at(0).pathLength(), straight.at(0).pathLength());
  }
  // At the plane but outside its finite bounds: prediction must not claim
  // arrival.
  const auto outside =
      plane->intersect(gctx, Vector3(11., 0., 20.), Vector3::UnitZ(),
                       Vector3(-0.005, 0., 0.), BoundaryTolerance::None());
  BOOST_CHECK(!outside.at(0).isValid());
  const auto inside =
      plane->intersect(gctx, Vector3(0., 0., 20.), Vector3::UnitZ(),
                       Vector3(0.005, 0., 0.), BoundaryTolerance::None());
  BOOST_CHECK(inside.at(0).status() == IntersectionStatus::onSurface);
}
BOOST_AUTO_TEST_CASE(CurvatureRequiresRaySeedAndFallsBackAtGrazingIncidence) {
  const auto plane = Surface::makeShared<PlaneSurface>(
      Transform3(Translation3(0., 0., 1e-8)),
      std::make_shared<RectangleBounds>(10., 10.));
  const Vector3 start = Vector3::Zero();
  const Vector3 curvature(0., 0., 0.001);
  const Vector3 grazing = Vector3(1., 0., 1e-7).normalized();
  // Use a tighter arrival tolerance so the start is not already on the plane.
  const double tolerance = 1e-12;
  const auto ray = plane
                       ->intersect(gctx, start, grazing,
                                   BoundaryTolerance::None(), tolerance)
                       .at(0);
  const auto fallback = plane
                            ->intersect(gctx, start, grazing, curvature,
                                        BoundaryTolerance::None(), tolerance)
                            .at(0);
  BOOST_REQUIRE(ray.isValid());
  CHECK_CLOSE_ABS(fallback.position(), ray.position(), 1e-12);
  CHECK_CLOSE_ABS(fallback.pathLength(), ray.pathLength(), 1e-12);
  const auto parallel =
      plane
          ->intersect(gctx, start, Vector3::UnitX(), curvature,
                      BoundaryTolerance::None(), tolerance)
          .at(0);
  BOOST_CHECK(!parallel.isValid());
}

BOOST_AUTO_TEST_CASE(CurvatureSurfaceDispatchAndDiscBounds) {
  const auto disc = Surface::makeShared<DiscSurface>(
      Transform3(Translation3(0., 0., 20.)),
      std::make_shared<RadialBounds>(0.5, 1.5, 0.2));
  const Surface& surface = *disc;
  const Vector3 origin = Vector3::Zero();
  const Vector3 direction = Vector3::UnitZ();
  BOOST_CHECK(
      !surface.intersect(gctx, origin, direction, BoundaryTolerance::None())
           .at(0)
           .isValid());
  const auto hit =
      surface
          .intersect(gctx, origin, direction, Vector3(0.005, 0., 0.),
                     BoundaryTolerance::None())
          .at(0);
  BOOST_REQUIRE(hit.isValid());
  CHECK_CLOSE_ABS(hit.position(), Vector3(1., 0., 20.), 1e-12);
  BOOST_CHECK(!disc->intersect(gctx, origin, direction, Vector3(0., 0.005, 0.),
                               BoundaryTolerance::None())
                   .at(0)
                   .isValid());
  BOOST_CHECK(!disc->intersect(gctx, origin, direction, Vector3(0.01, 0., 0.),
                               BoundaryTolerance::None())
                   .at(0)
                   .isValid());
  BOOST_CHECK(disc->intersect(gctx, origin, direction, Vector3(0.01, 0., 0.),
                              BoundaryTolerance::AbsoluteEuclidean(0.6))
                  .at(0)
                  .isValid());
  BOOST_CHECK(disc->intersect(gctx, origin, direction, Vector3(0., 0.005, 0.),
                              BoundaryTolerance::Infinite())
                  .at(0)
                  .isValid());

  // Unsupported shapes inherit the ray implementation; neither concrete nor
  // polymorphic calls require a curvature override.
  const auto point = Surface::makeShared<PointSurface>(Vector3(0., 0., 20.));
  const Surface& unsupported = *point;
  const auto ray = point->intersect(gctx, origin, direction).at(0);
  for (const auto& result :
       {unsupported.intersect(gctx, origin, direction, Vector3(0.005, 0., 0.)),
        point->intersect(gctx, origin, direction, Vector3(0.005, 0., 0.))}) {
    CHECK_CLOSE_ABS(result.at(0).position(), ray.position(), 1e-12);
    BOOST_CHECK_EQUAL(result.at(0).pathLength(), ray.pathLength());
  }
}

BOOST_AUTO_TEST_CASE(CurvatureCustomSurfaceKeepsRayImplementation) {
  const SurfaceStub custom;
  const Surface& surface = custom;
  const auto ray = surface.intersect(gctx, Vector3::Zero(), Vector3::UnitX());
  const auto result = surface.intersect(gctx, Vector3::Zero(), Vector3::UnitX(),
                                        Vector3(0., 0.001, 0.));
  BOOST_REQUIRE_EQUAL(result.size(), ray.size());
  CHECK_CLOSE_ABS(result.at(0).position(), ray.at(0).position(), 1e-12);
  BOOST_CHECK_EQUAL(result.at(0).pathLength(), ray.at(0).pathLength());
  BOOST_CHECK(result.at(1).status() == ray.at(1).status());
}

BOOST_AUTO_TEST_CASE(CurvatureCylinderBoundsAndArrival) {
  const auto cylinder =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20., 0.5);
  const Vector3 start = Vector3::Zero();
  const Vector3 direction = Vector3::UnitX();
  const Vector3 curvature(0., 0., 0.005);
  for (const auto& tolerance :
       {BoundaryTolerance::None(), BoundaryTolerance::AbsoluteEuclidean(0.6),
        BoundaryTolerance::Infinite()}) {
    const auto result =
        cylinder->intersect(gctx, start, direction, curvature, tolerance);
    BOOST_CHECK_EQUAL(result.size(), 2u);
    for (const auto& hit : result) {
      // Use the surface's existing tolerance semantics, including at the
      // azimuthal seam of the negative root.
      BOOST_CHECK_EQUAL(
          hit.isValid(),
          cylinder->isOnSurface(gctx, hit.position(), direction, tolerance));
      CHECK_CLOSE_ABS(hit.position().z(), 1., 1e-12);
      CHECK_CLOSE_ABS(hit.position().head<2>().norm(), 20., 1e-12);
    }
    BOOST_CHECK_EQUAL(result.at(1).isValid(), !tolerance.isNone());
    CHECK_CLOSE_ABS(result.at(0).pathLength(), -20., 1e-12);
    CHECK_CLOSE_ABS(result.at(1).pathLength(), 20., 1e-12);
  }
  const auto arrival =
      cylinder->intersect(gctx, Vector3(20., 0., 0.), direction, curvature,
                          BoundaryTolerance::None());
  BOOST_CHECK(arrival.at(1).status() == IntersectionStatus::onSurface);
  CHECK_CLOSE_ABS(arrival.at(1).position(), Vector3(20., 0., 0.), 1e-12);
  // A rejected correction must retain the original expanded ray-length bound.
  const Vector3 edge(0., 0., 0.5 + 0.5 * s_onSurfaceTolerance);
  const auto fallback = cylinder->intersect(
      gctx, edge, direction, Vector3(0., 1., 0.), BoundaryTolerance::None());
  BOOST_CHECK(fallback.at(0).isValid());
  BOOST_CHECK(fallback.at(1).isValid());
  // In contrast, accepted predictions check the actual requested bounds.
  const auto corrected = cylinder->intersect(
      gctx, edge, direction, Vector3(0., 0.001, 0.), BoundaryTolerance::None());
  BOOST_CHECK(!corrected.at(0).isValid());
  BOOST_CHECK(!corrected.at(1).isValid());
}

BOOST_AUTO_TEST_CASE(CurvatureInvalidInputAndTangentBranches) {
  const auto cylinder =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20., 10.);
  const Vector3 start(20., -5., 0.);
  const Vector3 direction = Vector3::UnitY();
  const auto ray =
      cylinder->intersect(gctx, start, direction, BoundaryTolerance::None());
  for (const Vector3& curvature :
       {Vector3::Zero().eval(), Vector3(0.001, 0., 0.),
        Vector3(std::numeric_limits<double>::quiet_NaN(), 0., 0.)}) {
    const auto result = cylinder->intersect(gctx, start, direction, curvature,
                                            BoundaryTolerance::None());
    BOOST_REQUIRE_EQUAL(result.size(), 2u);
    for (std::size_t i = 0; i < result.size(); ++i) {
      BOOST_CHECK(result.at(i).status() == ray.at(i).status());
      CHECK_CLOSE_ABS(result.at(i).position(), ray.at(i).position(), 1e-12);
      CHECK_CLOSE_ABS(result.at(i).pathLength(), ray.at(i).pathLength(), 1e-12);
    }
  }
}

}  // namespace ActsTests
