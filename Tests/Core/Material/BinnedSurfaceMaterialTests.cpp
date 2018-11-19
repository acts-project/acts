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
}
}
