// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/detail/CurvatureIntersection.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

using namespace Acts;
namespace ActsTests {
const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

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
  const auto crossing = detail::intersectWithCurvature(
      *sensor, gctx, start, direction, curvature, BoundaryTolerance::None());
  BOOST_REQUIRE(crossing.at(0).isValid());
  BOOST_CHECK_EQUAL(crossing.size(), 1u);
  CHECK_CLOSE_ABS(crossing.at(0).position(), Vector3(1., 0., 20.), 1e-12);
  CHECK_CLOSE_ABS(crossing.at(0).pathLength(), 20., 1e-12);
  // Analytic circle with initial tangent +z and curvature towards +x.
  const double exactX = (1. - std::sqrt(1. - std::pow(20. * 0.005, 2))) / 0.005;
  BOOST_CHECK_LT(std::abs(crossing.at(0).position().x() - exactX), 0.003);
  const auto reversed = detail::intersectWithCurvature(
      *sensor, gctx, start, -direction, curvature, BoundaryTolerance::None());
  CHECK_CLOSE_ABS(reversed.at(0).position(), crossing.at(0).position(), 1e-12);
  CHECK_CLOSE_ABS(reversed.at(0).pathLength(), -20., 1e-12);
  BOOST_CHECK(!detail::intersectWithCurvature(*sensor, gctx, start, direction,
                                              -curvature,
                                              BoundaryTolerance::None())
                   .at(0)
                   .isValid());

  Transform3 tilted = Transform3::Identity();
  tilted.translation() = Vector3(0., 0., 20.);
  tilted.linear() =
      AngleAxis3(std::atan(0.2), Vector3::UnitY()).toRotationMatrix();
  const auto plane = Surface::makeShared<PlaneSurface>(
      tilted, std::make_shared<RectangleBounds>(10., 10.));
  const auto hit =
      detail::intersectWithCurvature(*plane, gctx, start, direction, curvature,
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
  const auto hits = detail::intersectWithCurvature(
      *cylinder, gctx, Vector3::Zero(), direction, Vector3(0., 0.005, 0.),
      BoundaryTolerance::None());
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
  const auto movedHits = detail::intersectWithCurvature(
      *moved, gctx, transform * Vector3::Zero(), transform.linear() * direction,
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
    const auto fallback = detail::intersectWithCurvature(
        *plane, gctx, Vector3::Zero(), Vector3::UnitZ(), curvature,
        BoundaryTolerance::None());
    CHECK_CLOSE_ABS(fallback.at(0).position(), straight.at(0).position(),
                    1e-12);
    BOOST_CHECK_EQUAL(fallback.at(0).pathLength(), straight.at(0).pathLength());
  }
  // At the plane but outside its finite bounds: prediction must not claim
  // arrival.
  const auto outside = detail::intersectWithCurvature(
      *plane, gctx, Vector3(11., 0., 20.), Vector3::UnitZ(),
      Vector3(-0.005, 0., 0.), BoundaryTolerance::None());
  BOOST_CHECK(!outside.at(0).isValid());
  const auto inside = detail::intersectWithCurvature(
      *plane, gctx, Vector3(0., 0., 20.), Vector3::UnitZ(),
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
  const auto fallback =
      detail::intersectWithCurvature(*plane, gctx, start, grazing, curvature,
                                     BoundaryTolerance::None(), tolerance)
          .at(0);
  BOOST_REQUIRE(ray.isValid());
  CHECK_CLOSE_ABS(fallback.position(), ray.position(), 1e-12);
  CHECK_CLOSE_ABS(fallback.pathLength(), ray.pathLength(), 1e-12);
  const auto parallel = detail::intersectWithCurvature(
                            *plane, gctx, start, Vector3::UnitX(), curvature,
                            BoundaryTolerance::None(), tolerance)
                            .at(0);
  BOOST_CHECK(!parallel.isValid());
}

}  // namespace ActsTests
