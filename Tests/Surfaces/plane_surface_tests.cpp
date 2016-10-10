// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Plane Surface Tests
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace bdata = boost::unit_test::data;
namespace utf   = boost::unit_test;

namespace Acts {

namespace Test {
  /// Unit test for creating plane surface from normal vector and position
  ///
  /// Test for random plane positions and orientations.
  BOOST_TEST_DECORATOR(*utf::tolerance(1e-10))
  BOOST_DATA_TEST_CASE(plane_surface_normal1,
                       bdata::random(-10., 10.) ^ bdata::random(-10., 10.)
                           ^ bdata::random(-10., 10.)
                           ^ bdata::random(0., 2 * M_PI)
                           ^ bdata::random(0., M_PI)
                           ^ bdata::xrange(100),
                       x,
                       y,
                       z,
                       phi,
                       theta,
                       index)
  {
    const double   nx = sin(theta) * cos(phi);
    const double   ny = sin(theta) * sin(phi);
    const double   nz = cos(theta);
    const Vector3D pos(x, y, z);
    const Vector3D normal(nx, ny, nz);
    PlaneSurface   p(pos, normal);

    const auto& transform = p.transform();
    const auto& ex        = transform.matrix().col(0).head<3>();
    const auto& ey        = transform.matrix().col(1).head<3>();
    const auto& ez        = transform.matrix().col(2).head<3>();
    const auto& trans     = transform.matrix().col(3).head<3>();
    const auto& n         = p.normal();

    BOOST_TEST(n.dot(normal) == 1.);
    BOOST_TEST(n.dot(ex) == 0.);
    BOOST_TEST(n.dot(ey) == 0.);
    BOOST_TEST(n.dot(ez) == 1.);
    BOOST_TEST(ex.dot(ey) == 0.);
    BOOST_TEST(ex.norm() == 1.);
    BOOST_TEST(ey.norm() == 1.);
    BOOST_TEST(ez.norm() == 1.);
    BOOST_TEST(trans.dot(pos) == trans.norm() * pos.norm());
  }

  /// Unit test for creating plane surface from normal vector and position
  ///
  /// Test for random plane positions, but specific orientations to catch
  /// pathological situations.
  BOOST_TEST_DECORATOR(*utf::tolerance(1e-10))
  BOOST_DATA_TEST_CASE(plane_surface_normal2,
                       bdata::random(-10., 10.) ^ bdata::random(-10., 10.)
                           ^ bdata::random(-10., 10.)
                           ^ bdata::make(std::vector<float>({1, 0, 0}))
                           ^ bdata::make(std::vector<float>({0, 1, 0}))
                           ^ bdata::make(std::vector<float>({0, 0, 1}))
                           ^ bdata::xrange(3),
                       x,
                       y,
                       z,
                       nx,
                       ny,
                       nz,
                       index)
  {
    const Vector3D pos(x, y, z);
    const Vector3D normal(nx, ny, nz);
    PlaneSurface   p(pos, normal);

    const auto& transform = p.transform();
    const auto& ex        = transform.matrix().col(0).head<3>();
    const auto& ey        = transform.matrix().col(1).head<3>();
    const auto& ez        = transform.matrix().col(2).head<3>();
    const auto& trans     = transform.matrix().col(3).head<3>();
    const auto& n         = p.normal();

    BOOST_TEST(n.dot(normal) == 1.);
    BOOST_TEST(n.dot(ex) == 0.);
    BOOST_TEST(n.dot(ey) == 0.);
    BOOST_TEST(n.dot(ez) == 1.);
    BOOST_TEST(ex.dot(ey) == 0.);
    BOOST_TEST(ex.norm() == 1.);
    BOOST_TEST(ey.norm() == 1.);
    BOOST_TEST(ez.norm() == 1.);
    BOOST_TEST(trans.dot(pos) == trans.norm() * pos.norm());
  }
}  // end of namespace Test

}  // end of namespace Acts
