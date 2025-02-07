// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/Geometry/detail/PolyhedronIntersection.hpp>
#include <Acts/Surfaces/AnnulusBounds.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>

using namespace Acts;
using namespace Acts::detail;

BOOST_AUTO_TEST_CASE(coplanar_1) {
  std::array<Vector3, 3> t1 = {Vector3(0., 0., 0.), Vector3(1., 0., 0.),
                               Vector3(0., 1., 0.)};
  std::array<Vector3, 3> t2 = {Vector3(0., 0., 0.), Vector3(2., 0., 0.),
                               Vector3(0., 2., 0.)};
  BOOST_CHECK_EQUAL(triangleIntersection3D(t1, t2),
                    TriIntersectionResult::eTriCoplanar);
  BOOST_CHECK(triangleCoplanarIntersect(t1, t2) == true);
}

BOOST_AUTO_TEST_CASE(coplanar_2) {
  std::array<Vector3, 3> t1 = {Vector3(0., 0., 0.), Vector3(1., 0., 0.),
                               Vector3(0., 1., 0.)};
  // shift by five in x
  std::array<Vector3, 3> t2 = {Vector3(5., 0., 0.), Vector3(6., 0., 0.),
                               Vector3(5., 1., 0.)};
  BOOST_CHECK_EQUAL(triangleIntersection3D(t1, t2),
                    TriIntersectionResult::eTriCoplanar);
  BOOST_CHECK(triangleCoplanarIntersect(t1, t2) == false);
}

BOOST_AUTO_TEST_CASE(shifted_z) {
  using namespace Acts::detail;
  std::array<Vector3, 3> t1 = {Vector3(0., 0., 1.), Vector3(1., 0., 1.),
                               Vector3(0., 1., 1.)};
  std::array<Vector3, 3> t2 = {Vector3(0., 0., 0.), Vector3(1., 0., 0.),
                               Vector3(1., 1., 0.)};
  BOOST_CHECK_EQUAL(triangleIntersection3D(t1, t2),
                    TriIntersectionResult::eTriIntersectFalse);
}

BOOST_AUTO_TEST_CASE(shifted_x) {
  using namespace Acts::detail;
  std::array<Vector3, 3> t1 = {Vector3(1., 0., 0.), Vector3(1., 1., 0.),
                               Vector3(1., 0., 1.)};
  std::array<Vector3, 3> t2 = {Vector3(0., 0., 0.), Vector3(0., 1., 0.),
                               Vector3(0., 0., 1.)};
  BOOST_CHECK_EQUAL(triangleIntersection3D(t1, t2),
                    TriIntersectionResult::eTriIntersectFalse);
}

BOOST_AUTO_TEST_CASE(cross_intersect_1) {
  using namespace Acts::detail;
  std::array<Vector3, 3> t1 = {Vector3(0., 0., 0.), Vector3(1., 0., 0.),
                               Vector3(0., 1., 0.)};
  std::array<Vector3, 3> t2 = {Vector3(0., 0.5, 0.), Vector3(1., 0.5, 0.),
                               Vector3(0., 0.5, 1.)};
  BOOST_CHECK_EQUAL(triangleIntersection3D(t1, t2),
                    TriIntersectionResult::eTriIntersectTrue);
}

BOOST_AUTO_TEST_CASE(cross_intersect_2) {
  using namespace Acts::detail;
  std::array<Vector3, 3> t1 = {Vector3(0., 0., 0.), Vector3(1., 0., 0.),
                               Vector3(0., 1., 0.)};
  std::array<Vector3, 3> t2 = {Vector3(0., 2., 0.), Vector3(1., 2., 0.),
                               Vector3(0., 2., 1.)};
  BOOST_CHECK_EQUAL(triangleIntersection3D(t1, t2),
                    TriIntersectionResult::eTriIntersectFalse);
}

template <typename surface_t, typename bounds_t>
bool checkSurfaces(Transform3 trafo1, Transform3 trafo2,
                   std::shared_ptr<bounds_t> bounds, std::size_t nPoints = 4) {
  auto surface1 = Surface::makeShared<surface_t>(trafo1, bounds);
  auto surface2 = Surface::makeShared<surface_t>(trafo2, bounds);

  auto poly1 = surface1->polyhedronRepresentation({}, nPoints);
  // BOOST_CHECK(poly1.vertices.size() == nPoints);
  // BOOST_CHECK(poly1.triangularMesh.size() >= 2);
  auto poly2 = surface2->polyhedronRepresentation({}, nPoints);
  // BOOST_CHECK(poly2.vertices.size() == nPoints);
  // BOOST_CHECK(poly2.triangularMesh.size() >= 2);

  return polyhedronIntersection(poly1, poly2);
}

BOOST_AUTO_TEST_CASE(polyhedron_intersection_test_1) {
  auto bounds = std::make_shared<RectangleBounds>(1.0, 1.0);
  auto trafo1 = Transform3::Identity();
  auto trafo2 = Transform3::Identity();
  BOOST_CHECK((checkSurfaces<PlaneSurface, RectangleBounds>(
                   trafo1, trafo2, bounds, 4) == true));
}

BOOST_AUTO_TEST_CASE(polyhedron_intersection_test_2) {
  auto bounds = std::make_shared<RectangleBounds>(1.0, 1.0);
  auto trafo1 = Transform3::Identity();
  auto trafo2 = Transform3::Identity();
  trafo2.translate(Vector3(10., 0., 0.));
  BOOST_CHECK((checkSurfaces<PlaneSurface, RectangleBounds>(
                   trafo1, trafo2, bounds, 4) == false));
}

const double minRadius = 7.2;
const double maxRadius = 12.0;
const double minPhi = 0.74195;
const double maxPhi = 1.33970;
const Vector2 offset(-2., 2.);

BOOST_AUTO_TEST_CASE(polyhedron_intersection_test_annulus_1) {
  std::shared_ptr<DiscBounds> bounds = std::make_shared<AnnulusBounds>(
      minRadius, maxRadius, minPhi, maxPhi, offset);
  auto trafo1 = Transform3::Identity();
  auto trafo2 = Transform3::Identity();
  BOOST_CHECK((checkSurfaces<DiscSurface, DiscBounds>(trafo1, trafo2, bounds,
                                                      36) == true));
}

BOOST_AUTO_TEST_CASE(polyhedron_intersection_test_annulus_2) {
  std::shared_ptr<DiscBounds> bounds = std::make_shared<AnnulusBounds>(
      minRadius, maxRadius, minPhi, maxPhi, offset);
  auto trafo1 = Transform3::Identity();
  auto trafo2 = Transform3::Identity();
  trafo2.translate(Vector3(10.0, 0., 0.));
  BOOST_CHECK((checkSurfaces<DiscSurface, DiscBounds>(trafo1, trafo2, bounds,
                                                      36) == false));
}

BOOST_AUTO_TEST_CASE(pathological_1) {
  auto t1 = std::array{Vector3(31.6775, -9.31483, -243.45),
                       Vector3(36.3225, 9.31483, -243.45),
                       Vector3(36.3225, 9.31483, -223.45)};
  auto t2 = std::array{Vector3(31.6775, -9.31483, -223.15),
                       Vector3(36.3225, 9.31483, -203.15),
                       Vector3(31.6775, -9.31483, -203.15)};

  BOOST_CHECK(triangleCoplanarIntersect(t1, t2, 1e-3) == false);
}
