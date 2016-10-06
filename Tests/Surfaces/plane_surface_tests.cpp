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
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {
  /// Unit test for creating plane surface from normal vector and position
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

    BOOST_TEST(p.normal().dot(normal) == 1., tt::tolerance(1e-10));
  }

  /// Unit test for creating plane surface from normal vector and position
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

    BOOST_TEST(p.normal().dot(normal) == 1., tt::tolerance(1e-10));
  }
}  // end of namespace Test

}  // end of namespace Acts
