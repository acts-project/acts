// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(SurfacesSuite)
/// Unit test for creating compliant/non-compliant InfiniteBounds object
BOOST_AUTO_TEST_CASE(InfiniteBoundsConstruction) {
  InfiniteBounds u;
  BOOST_CHECK_EQUAL(u.type(), SurfaceBounds::Boundless);
  // InfiniteBounds s(1);  // would act as std::size_t cast to InfiniteBounds
  // InfiniteBounds t(s);
  InfiniteBounds v(u);  // implicit
  BOOST_CHECK_EQUAL(v.type(), SurfaceBounds::Boundless);
}
/// Unit tests for InfiniteBounds properties
BOOST_AUTO_TEST_CASE(InfiniteBoundsProperties) {
  InfiniteBounds infiniteBoundsObject;
  /// Test for type()
  BOOST_CHECK_EQUAL(infiniteBoundsObject.type(), SurfaceBounds::Boundless);

  /// Test for inside()
  const Vector2 anyVector{0., 1.};
  const BoundaryTolerance anyTolerance = BoundaryTolerance::None();
  BOOST_CHECK(infiniteBoundsObject.inside(anyVector, anyTolerance));

  /// Test for dump
  boost::test_tools::output_test_stream dumpOutput;
  infiniteBoundsObject.toStream(dumpOutput);
  BOOST_CHECK(
      dumpOutput.is_equal("Acts::InfiniteBounds ... boundless surface\n"));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
