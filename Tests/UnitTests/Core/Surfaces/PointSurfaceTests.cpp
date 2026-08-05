// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/PointSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <memory>
#include <numbers>
#include <tuple>
#include <vector>

using namespace Acts;

namespace ActsTests {

// Create a test context
GeometryContext tgContext = GeometryContext::dangerouslyDefaultConstruct();

BOOST_AUTO_TEST_SUITE(SurfacesSuite)

/// Construction, type, name, bounds
BOOST_AUTO_TEST_CASE(PointSurfaceConstruction) {
  Vector3 center{1., 2., 3.};

  // unbounded
  auto unbounded = Surface::makeShared<PointSurface>(center);
  BOOST_CHECK_EQUAL(unbounded->type(), Surface::Point);
  BOOST_CHECK_EQUAL(unbounded->name(), "Acts::PointSurface");
  CHECK_CLOSE_ABS(unbounded->center(tgContext), center, 1e-12);
  // unbounded surface reports the shared boundless bounds
  BOOST_CHECK_EQUAL(unbounded->bounds().type(), SurfaceBounds::eBoundless);
  BOOST_CHECK(unbounded->boundsPtr() == nullptr);

  // bounded
  auto bounded = Surface::makeShared<PointSurface>(center, 5.);
  BOOST_CHECK_EQUAL(bounded->type(), Surface::Point);
  BOOST_CHECK_EQUAL(bounded->bounds().type(), SurfaceBounds::ePoint);

  // copy
  auto copied = Surface::makeShared<PointSurface>(*bounded);
  BOOST_CHECK_EQUAL(copied->type(), Surface::Point);
  BOOST_CHECK(*copied == *bounded);
}

/// The measurement frame normal equals the (normalized) direction
BOOST_AUTO_TEST_CASE(PointSurfaceNormalIsDirection) {
  auto surface = Surface::makeShared<PointSurface>(Vector3(0., 0., 0.));
  Vector3 dir = Vector3(1., 2., 3.).normalized();
  CHECK_CLOSE_ABS(surface->normal(tgContext, Vector3::Zero(), dir), dir, 1e-12);

  auto rframe = surface->referenceFrame(tgContext, Vector3::Zero(), dir);
  CHECK_CLOSE_ABS(rframe.col(2), dir, 1e-12);
  // orthonormal frame
  CHECK_CLOSE_ABS(rframe.col(0).dot(rframe.col(1)), 0., 1e-12);
  CHECK_CLOSE_ABS(rframe.col(0).norm(), 1., 1e-12);
  CHECK_CLOSE_ABS(rframe.col(1).norm(), 1., 1e-12);
}

/// Round-trip: intersect -> globalToLocal -> localToGlobal
BOOST_AUTO_TEST_CASE(PointSurfaceRoundTrip) {
  Vector3 center{2., -1., 0.5};
  auto surface = Surface::makeShared<PointSurface>(center);

  auto roundTrip = [&](const Vector3& pos, const Vector3& dir) {
    Intersection3D intersection =
        surface->intersect(tgContext, pos, dir).closest();
    Vector3 global = intersection.position();
    Vector2 local = *surface->globalToLocal(tgContext, global, dir);
    Vector3 global2 = surface->localToGlobal(tgContext, local, dir);
    return std::make_tuple(global, local, global2);
  };

  {
    Vector3 pos = {-0.5, 0.3, 1.2};
    Vector3 dir = Vector3(-0.2, -0.7, 0.4).normalized();
    auto [global, local, global2] = roundTrip(pos, dir);
    CHECK_CLOSE_ABS(global, global2, 1e-10);
  }
  {
    Vector3 pos = {5., 5., 5.};
    Vector3 dir = Vector3(1., -1., 0.2).normalized();
    auto [global, local, global2] = roundTrip(pos, dir);
    CHECK_CLOSE_ABS(global, global2, 1e-10);
  }
}

/// Round-trip stability including large eta and direction parallel to Z
BOOST_AUTO_TEST_CASE(PointSurfaceRoundTripEtaStability) {
  Vector3 center = {3., 0., 0.};
  auto surface = Surface::makeShared<PointSurface>(center);

  const std::vector<double> etas = {0, 1, 2, 3, 5, 10};

  for (double eta : etas) {
    Vector3 dir = makeDirectionFromPhiEta(std::numbers::pi / 4., eta);
    Vector3 pos = center + Vector3(0.7, -0.4, 0.9);

    Intersection3D intersection =
        surface->intersect(tgContext, pos, dir).closest();
    Vector3 global = intersection.position();
    Vector2 local = *surface->globalToLocal(tgContext, global, dir);
    Vector3 global2 = surface->localToGlobal(tgContext, local, dir);

    CHECK_CLOSE_ABS(global, global2, 1e-9);
    // the PCA residual must be perpendicular to the direction
    CHECK_CLOSE_ABS((global - center).dot(dir), 0., 1e-9);
  }

  // direction (anti-)parallel to Z (non-standard frame branch)
  for (Vector3 dir : {Vector3(0., 0., 1.), Vector3(0., 0., -1.)}) {
    Vector3 pos = center + Vector3(1., 2., 4.);
    Intersection3D intersection =
        surface->intersect(tgContext, pos, dir).closest();
    Vector3 global = intersection.position();
    Vector2 local = *surface->globalToLocal(tgContext, global, dir);
    Vector3 global2 = surface->localToGlobal(tgContext, local, dir);
    CHECK_CLOSE_ABS(global, global2, 1e-9);
    CHECK_CLOSE_ABS((global - center).dot(dir), 0., 1e-9);
  }
}

/// The intersection is the point of closest approach; no degeneracy
BOOST_AUTO_TEST_CASE(PointSurfaceIntersection) {
  Vector3 center{1., 1., 1.};
  auto surface = Surface::makeShared<PointSurface>(center);

  Vector3 pos{0., 0., 0.};
  // Try several directions, including axis-aligned ones (would be degenerate
  // for a LineSurface but are fine here).
  const std::vector<Vector3> dirs = {
      Vector3(1., 0., 0.), Vector3(0., 1., 0.), Vector3(0., 0., 1.),
      Vector3(1., 1., 1.).normalized(), Vector3(-1., 2., -0.5).normalized()};

  for (const Vector3& dir : dirs) {
    Intersection3D intersection =
        surface->intersect(tgContext, pos, dir).closest();
    BOOST_CHECK(intersection.isValid());

    double u = (center - pos).dot(dir);
    Vector3 expected = pos + u * dir;
    CHECK_CLOSE_ABS(intersection.position(), expected, 1e-10);
    CHECK_CLOSE_ABS(intersection.pathLength(), u, 1e-10);
    // residual perpendicular to the direction
    CHECK_CLOSE_ABS((intersection.position() - center).dot(dir), 0., 1e-10);
  }

  // starting exactly at the PCA reports onSurface
  Vector3 dir = Vector3(1., 0., 0.);
  Vector3 pca = pos + (center - pos).dot(dir) * dir;
  Intersection3D onSurf = surface->intersect(tgContext, pca, dir).closest();
  BOOST_CHECK_EQUAL(onSurf.status(), IntersectionStatus::onSurface);
}

/// Bounds are respected in the intersection boundary check
BOOST_AUTO_TEST_CASE(PointSurfaceIntersectionBounds) {
  Vector3 center{0., 0., 0.};
  // maximum distance 1
  auto bounded = Surface::makeShared<PointSurface>(center, 1.);
  auto unbounded = Surface::makeShared<PointSurface>(center);

  // A track passing at transverse distance ~2 from the point
  Vector3 pos{2., 0., 0.};
  Vector3 dir = Vector3(0., 0., 1.);

  auto within = BoundaryTolerance::None();

  Intersection3D boundedHit =
      bounded->intersect(tgContext, pos, dir, within).closest();
  BOOST_CHECK_EQUAL(boundedHit.status(), IntersectionStatus::unreachable);

  Intersection3D unboundedHit =
      unbounded->intersect(tgContext, pos, dir, within).closest();
  BOOST_CHECK(unboundedHit.isValid());

  // A track passing close to the point is inside the bounds
  Vector3 posClose{0.3, 0., 0.};
  Intersection3D closeHit =
      bounded->intersect(tgContext, posClose, dir, within).closest();
  BOOST_CHECK(closeHit.isValid());
}

// Build a free vector from a bound vector consistent with the surface frame
FreeVector boundToFree(const PointSurface& surface, const BoundVector& b) {
  Vector3 dir = makeDirectionFromPhiTheta(b[eBoundPhi], b[eBoundTheta]);
  Vector3 pos = surface.localToGlobal(
      tgContext, Vector2(b[eBoundLoc0], b[eBoundLoc1]), dir);
  FreeVector f = FreeVector::Zero();
  f.segment<3>(eFreePos0) = pos;
  f[eFreeTime] = b[eBoundTime];
  f.segment<3>(eFreeDir0) = dir;
  f[eFreeQOverP] = b[eBoundQOverP];
  return f;
}

/// boundToFreeJacobian validated against finite differences at nonzero local
BOOST_AUTO_TEST_CASE(PointSurfaceBoundToFreeJacobian) {
  auto surface = Surface::makeShared<PointSurface>(Vector3(1., -2., 0.5));

  // pick a non-forward direction (standard frame branch) and nonzero local so
  // the frame-rotation correction is exercised
  double phi = 0.6;
  double theta = std::numbers::pi / 3.;
  BoundVector b;
  b << 0.7, -0.4, phi, theta, 0.5, 1.5;

  FreeVector f0 = boundToFree(*surface, b);
  Vector3 pos0 = f0.segment<3>(eFreePos0);
  Vector3 dir0 = f0.segment<3>(eFreeDir0);

  BoundToFreeMatrix jac = surface->boundToFreeJacobian(tgContext, pos0, dir0);

  double h = 1e-6;
  for (int i = 0; i < static_cast<int>(eBoundSize); ++i) {
    BoundVector bp = b;
    BoundVector bm = b;
    bp[i] += h;
    bm[i] -= h;
    FreeVector numeric =
        (boundToFree(*surface, bp) - boundToFree(*surface, bm)) / (2 * h);
    for (int j = 0; j < static_cast<int>(eFreeSize); ++j) {
      BOOST_CHECK_SMALL(numeric[j] - jac(j, i), 1e-6);
    }
  }
}

/// freeToPathDerivative validated against finite differences
BOOST_AUTO_TEST_CASE(PointSurfaceFreeToPathDerivative) {
  Vector3 center{1., -2., 0.5};
  auto surface = Surface::makeShared<PointSurface>(center);

  Vector3 dir = Vector3(0.3, -0.5, 0.8).normalized();
  // start on the surface (PCA) so the assert holds
  Vector3 pos = center + Vector3(0.4, 0.1, -0.2);
  // project pos onto the measurement plane through the center
  pos -= (pos - center).dot(dir) * dir;

  FreeToPathMatrix analytic =
      surface->freeToPathDerivative(tgContext, pos, dir);

  auto pathU = [&](const Vector3& p, const Vector3& d) {
    return (center - p).dot(d);
  };

  double h = 1e-7;
  for (int k = 0; k < 3; ++k) {
    Vector3 pp = pos, pm = pos;
    pp[k] += h;
    pm[k] -= h;
    double numeric = (pathU(pp, dir) - pathU(pm, dir)) / (2 * h);
    BOOST_CHECK_SMALL(numeric - analytic(0, eFreePos0 + k), 1e-6);
  }
  for (int k = 0; k < 3; ++k) {
    Vector3 dp = dir, dm = dir;
    dp[k] += h;
    dm[k] -= h;
    double numeric = (pathU(pos, dp) - pathU(pos, dm)) / (2 * h);
    BOOST_CHECK_SMALL(numeric - analytic(0, eFreeDir0 + k), 1e-6);
  }
}

/// alignmentToPathDerivative: only the center contributes (= direction)
BOOST_AUTO_TEST_CASE(PointSurfaceAlignmentToPathDerivative) {
  Vector3 center{1., -2., 0.5};
  auto surface = Surface::makeShared<PointSurface>(center);
  Vector3 dir = Vector3(0.3, -0.5, 0.8).normalized();
  Vector3 pos = center + Vector3(0.4, 0.1, -0.2);
  pos -= (pos - center).dot(dir) * dir;

  AlignmentToPathMatrix analytic =
      surface->alignmentToPathDerivative(tgContext, pos, dir);

  // center block equals the direction
  CHECK_CLOSE_ABS(analytic.segment<3>(eAlignmentCenter0).transpose(), dir,
                  1e-10);
  // rotation block is zero
  CHECK_CLOSE_ABS(analytic.segment<3>(eAlignmentRotation0).norm(), 0., 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
