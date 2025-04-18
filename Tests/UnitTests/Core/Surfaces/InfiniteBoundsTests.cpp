// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)
/// Unit test for creating compliant/non-compliant InfiniteBounds object
BOOST_AUTO_TEST_CASE(InfiniteBoundsConstruction) {
  InfiniteBounds u;
  BOOST_CHECK_EQUAL(u.type(), SurfaceBounds::eBoundless);
  // InfiniteBounds s(1);  // would act as std::size_t cast to InfiniteBounds
  // InfiniteBounds t(s);
  InfiniteBounds v(u);  // implicit
  BOOST_CHECK_EQUAL(v.type(), SurfaceBounds::eBoundless);
}
/// Unit tests for InfiniteBounds properties
BOOST_AUTO_TEST_CASE(InfiniteBoundsProperties) {
  InfiniteBounds infiniteBoundsObject;
  /// test for type()
  BOOST_CHECK_EQUAL(infiniteBoundsObject.type(), SurfaceBounds::eBoundless);

  /// test for inside()
  const Vector2 anyVector{0., 1.};
  const BoundaryCheck anyBoundaryCheck(true);
  BOOST_CHECK(infiniteBoundsObject.inside(anyVector, anyBoundaryCheck));

  /// test for dump
  boost::test_tools::output_test_stream dumpOuput;
  infiniteBoundsObject.toStream(dumpOuput);
  BOOST_CHECK(
      dumpOuput.is_equal("Acts::InfiniteBounds ... boundless surface\n"));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
