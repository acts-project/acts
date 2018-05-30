// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Line Bounds Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit test for creating compliant/non-compliant LineBounds object
  BOOST_AUTO_TEST_CASE(LineBoundsConstruction)
  {
    /// test default construction
    LineBounds defaultConstructedLineBounds;  // implicit
    BOOST_TEST(defaultConstructedLineBounds.type() == SurfaceBounds::Line,
               "Test trivial method on default constructed object");
    /// test LineBounds(double, double)
    double radius(0.5), halfz(10.);
    BOOST_TEST(LineBounds(radius, halfz).type() == SurfaceBounds::Line);
    //
    LineBounds s(1);  // would act as size_t cast to LineBounds
    /// test copy construction;
    LineBounds copyConstructedLineBounds(
        defaultConstructedLineBounds);  // implicit
    BOOST_TEST(copyConstructedLineBounds.type() == SurfaceBounds::Line,
               "Test trivial method on copy constructed object");
  }

  /// Unit tests for LineBounds properties
  BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(LineBoundsProperties, 1)
  BOOST_AUTO_TEST_CASE(LineBoundsProperties)
  {
    // LineBounds object of radius 0.5 and halfz 20
    double     nominalRadius{0.5};
    double     nominalHalfLength{20.};
    LineBounds lineBoundsObject(nominalRadius, nominalHalfLength);

    /// test for clone
    auto pLineBoundsClone = lineBoundsObject.clone();
    BOOST_CHECK(pLineBoundsClone != nullptr);
    delete pLineBoundsClone;

    /// test for type()
    BOOST_TEST(lineBoundsObject.type() == SurfaceBounds::Line);

    /// test for inside()
    const Vector2D      origin{0., 0.};
    const Vector2D      atRadius{0.5, 10.};
    const Vector2D      beyondEnd{0.0, 30.0};
    const Vector2D      unitZ{0.0, 1.0};
    const Vector2D      unitR{1.0, 0.0};
    const BoundaryCheck trueBoundaryCheckWithTolerance(true, true, 0.1, 0.1);
    BOOST_TEST(lineBoundsObject.inside(atRadius, trueBoundaryCheckWithTolerance)
               == true);

    /// test for distanceToBoundary
    BOOST_TEST(lineBoundsObject.distanceToBoundary(unitR) == 1.);  // why?

    /// test for r()
    BOOST_TEST(lineBoundsObject.r() == nominalRadius);

    /// test for halflengthZ (NOTE: Naming violation)
    BOOST_TEST(lineBoundsObject.halflengthZ() == nominalHalfLength);

    /// test for dump
    boost::test_tools::output_test_stream dumpOuput;
    lineBoundsObject.dump(dumpOuput);
    BOOST_TEST(dumpOuput.is_equal(
        "Acts::LineBounds: (radius, halflengthInZ) = (0.5000000, 20.0000000)"));
  }
  /// Unit test for testing LineBounds assignment
  BOOST_AUTO_TEST_CASE(LineBoundsAssignment)
  {
    double     nominalRadius{0.5};
    double     nominalHalfLength{20.};
    LineBounds lineBoundsObject(nominalRadius, nominalHalfLength);
    LineBounds assignedLineBounds;
    assignedLineBounds = lineBoundsObject;
    BOOST_TEST(assignedLineBounds.r() == lineBoundsObject.r());
    BOOST_TEST(assignedLineBounds.halflengthZ()
               == lineBoundsObject.halflengthZ());
  }
  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
