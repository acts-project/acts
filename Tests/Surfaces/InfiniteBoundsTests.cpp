// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Infinite Bounds Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit test for creating compliant/non-compliant InfiniteBounds object
  BOOST_AUTO_TEST_CASE(InfiniteBoundsConstruction)
  {
    InfiniteBounds u;
    BOOST_TEST(u.type() == SurfaceBounds::Boundless,
               "Test trivial method on default constructed object");
    // InfiniteBounds s(1);  // would act as size_t cast to InfiniteBounds
    // InfiniteBounds t(s);
    InfiniteBounds v(u);  // implicit
    BOOST_TEST(v.type() == SurfaceBounds::Boundless,
               "Test trivial method on copy constructed object");
  }
  /// Unit tests for InfiniteBounds properties
  BOOST_AUTO_TEST_CASE(InfiniteBoundsProperties)
  {
    InfiniteBounds infiniteBoundsObject;
    /// test for type()
    BOOST_TEST(infiniteBoundsObject.type() == SurfaceBounds::Boundless);

    /// test for inside()
    const Vector2D      anyVector{0., 1.};
    const BoundaryCheck anyBoundaryCheck(true);
    BOOST_TEST(infiniteBoundsObject.inside(anyVector, anyBoundaryCheck)
               == true);

    /// test for distanceToBoundary
    BOOST_TEST_MESSAGE("Perhaps the following should be inf?");
    BOOST_TEST(infiniteBoundsObject.distanceToBoundary(anyVector) == 0.);

    /// test for clone
    auto pInfiniteBoundsClone = infiniteBoundsObject.clone();
    BOOST_CHECK(pInfiniteBoundsClone != nullptr);
    delete pInfiniteBoundsClone;

    /// test for dump
    boost::test_tools::output_test_stream dumpOuput;
    infiniteBoundsObject.dump(dumpOuput);
    BOOST_TEST(
        dumpOuput.is_equal("Acts::InfiniteBounds ... boundless surface\n"));
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
