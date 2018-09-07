// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE SurfaceMaterial Tests
#include <boost/test/included/unit_test.hpp>
#include <climits>
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"

namespace Acts {

namespace Test {

  /// Test the constructors
  BOOST_AUTO_TEST_CASE(MaterialProperties_construction_test)
  {
    // constructor only from argumnets
    MaterialProperties a(1., 2., 3., 4., 5., 6.);
    /// constructor with material
    MaterialProperties b(Material(1., 2., 3., 4., 5.), 6.);
    /// Check if they are equal
    BOOST_CHECK_EQUAL(a, b);

    /// Check the move construction
    MaterialProperties bMoved(std::move(b));
    /// Check if they are equal
    BOOST_CHECK_EQUAL(a, bMoved);

    /// Check the move assignment
    MaterialProperties bMovedAssigned = std::move(bMoved);
    /// Check if they are equal
    BOOST_CHECK_EQUAL(a, bMovedAssigned);
  }

  BOOST_AUTO_TEST_CASE(MaterialProperties_adding_test)
  {

    MaterialProperties a(1., 2., 3., 4., 5., 6.);
    MaterialProperties c(1., 2., 3., 4., 5., 0.);
    MaterialProperties d(0., 0., 0., 0., 0., 1.);

    // The average of twice the same material a should be a again
    a.average(a);
    BOOST_CHECK_EQUAL(a, a);
    // Adding material with 0 thickness should not change anything
    a.average(c);
    BOOST_CHECK_EQUAL(a, a);
    // Adding material with not material paramters (vacuum) should not change
    // anything
    a.average(d);
    BOOST_CHECK_EQUAL(a, a);
    // Adding material to no material paramters (vacuum) should give the added
    // material
    d.average(a);
    BOOST_CHECK_EQUAL(d, a);
  }
}
}
