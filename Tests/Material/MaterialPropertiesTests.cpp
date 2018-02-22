// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Material Tests
#include <boost/test/included/unit_test.hpp>
#include <climits>
#include "ACTS/Material/MaterialProperties.hpp"

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_CASE(AddMaterialProperties_test)
  {

    MaterialProperties a(1., 2., 3., 4., 5., 6.);
    MaterialProperties c(1., 2., 3., 4., 5., 0.);
    MaterialProperties d(0., 0., 0., 0., 0., 1.);

    // The average of twice the same material a should be a again
    a.add(a);
    BOOST_CHECK_EQUAL(a, a);
    // Adding material with 0 thickness should not change anything
    a.add(c);
    BOOST_CHECK_EQUAL(a, a);
    // Adding material with not material paramters (vacuum) should not change
    // anything
    a.add(d);
    BOOST_CHECK_EQUAL(a, a);
    // Adding material to no material paramters (vacuum) should give the added
    // material
    d.add(a);
    BOOST_CHECK_EQUAL(d, a);
  }
}
}
