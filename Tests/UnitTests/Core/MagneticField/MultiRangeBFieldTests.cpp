// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldError.hpp"
#include "Acts/MagneticField/MultiRangeBField.hpp"
#include "Acts/Utilities/Result.hpp"

using namespace Acts;

namespace ActsTests {

MagneticFieldContext mfContext = MagneticFieldContext();

BOOST_AUTO_TEST_SUITE(MagneticFieldSuite)

BOOST_AUTO_TEST_CASE(TestMultiRangeBField) {
  std::vector<std::pair<RangeXD<3, double>, Vector3>> inputs;

  inputs.emplace_back(RangeXD<3, double>{{0., 0., 0.}, {3., 3., 3.}},
                      Vector3{0., 0., 2.});
  inputs.emplace_back(RangeXD<3, double>{{1., 1., 1.}, {2., 2., 10.}},
                      Vector3{2., 0., 0.});

  const MultiRangeBField bfield(std::move(inputs));

  auto bcache = bfield.makeCache(mfContext);

  // Test a point inside the first volume.
  {
    Result<Vector3> r = bfield.getField({0.5, 0.5, 0.5}, bcache);
    BOOST_CHECK(r.ok());
    BOOST_CHECK_CLOSE((*r)[0], 0., 0.05);
    BOOST_CHECK_CLOSE((*r)[1], 0., 0.05);
    BOOST_CHECK_CLOSE((*r)[2], 2., 0.05);
  }

  // Test a point inside the second volume, not overlapping the first.
  {
    Result<Vector3> r = bfield.getField({1.5, 1.5, 5.}, bcache);
    BOOST_CHECK(r.ok());
    BOOST_CHECK_CLOSE((*r)[0], 2., 0.05);
    BOOST_CHECK_CLOSE((*r)[1], 0., 0.05);
    BOOST_CHECK_CLOSE((*r)[2], 0., 0.05);
  }

  // Test a point inside the first volume again.
  {
    Result<Vector3> r = bfield.getField({2.5, 2.5, 2.5}, bcache);
    BOOST_CHECK(r.ok());
    BOOST_CHECK_CLOSE((*r)[0], 0., 0.05);
    BOOST_CHECK_CLOSE((*r)[1], 0., 0.05);
    BOOST_CHECK_CLOSE((*r)[2], 2., 0.05);
  }

  // Test a point inside the second volume in the overlap region.
  {
    Result<Vector3> r = bfield.getField({1.5, 1.5, 1.5}, bcache);
    BOOST_CHECK(r.ok());
    BOOST_CHECK_CLOSE((*r)[0], 2., 0.05);
    BOOST_CHECK_CLOSE((*r)[1], 0., 0.05);
    BOOST_CHECK_CLOSE((*r)[2], 0., 0.05);
  }

  // Test a point outside all volumes.
  {
    Result<Vector3> r = bfield.getField({-1., -1., -1}, bcache);
    BOOST_CHECK(!r.ok());
  }
}

BOOST_AUTO_TEST_CASE(TestMultiRangeBFieldGradient) {
  std::vector<std::pair<RangeXD<3, double>, Vector3>> inputs;
  inputs.emplace_back(RangeXD<3, double>{{0., 0., 0.}, {3., 3., 3.}},
                      Vector3{0., 0., 2.});
  inputs.emplace_back(RangeXD<3, double>{{1., 1., 1.}, {2., 2., 2.}},
                      Vector3{2., 0., 0.});

  const MultiRangeBField bfield(std::move(inputs));
  auto bcache = bfield.makeCache(mfContext);

  BOOST_CHECK(bfield.providesFieldGradient());

  auto inside = bfield.getFieldAndGradient({0.5, 0.5, 0.5}, bcache);
  BOOST_REQUIRE(inside.ok());
  BOOST_CHECK_EQUAL(inside->field, Vector3(0., 0., 2.));
  BOOST_CHECK_EQUAL(inside->gradient, SquareMatrix3::Zero());

  // The second region overlaps the first one and has the higher priority
  auto overlap = bfield.getFieldAndGradient({1.5, 1.5, 1.5}, bcache);
  BOOST_REQUIRE(overlap.ok());
  BOOST_CHECK_EQUAL(overlap->field, Vector3(2., 0., 0.));
  BOOST_CHECK_EQUAL(overlap->gradient, SquareMatrix3::Zero());

  // The field must agree with getField at each position
  for (const Vector3& position :
       {Vector3(0.5, 0.5, 0.5), Vector3(1.5, 1.5, 1.5), Vector3(2.5, 1.5, 1.5),
        Vector3(1.5, 0.5, 1.5)}) {
    auto field = bfield.getField(position, bcache);
    auto fieldAndGradient = bfield.getFieldAndGradient(position, bcache);
    BOOST_REQUIRE(field.ok());
    BOOST_REQUIRE(fieldAndGradient.ok());
    BOOST_CHECK_EQUAL(fieldAndGradient->field, *field);
  }

  auto outside = bfield.getFieldAndGradient({5., 5., 5.}, bcache);
  BOOST_REQUIRE(!outside.ok());
  BOOST_CHECK(outside.error() == MagneticFieldError::OutOfBounds);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
