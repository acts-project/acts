// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE StrawSurface Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//
#include <limits>
#include "ACTS/Layers/PlaneLayer.hpp"
#include "ACTS/Material/HomogeneousSurfaceMaterial.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/StrawSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "DetectorElementStub.hpp"

using boost::test_tools::output_test_stream;
namespace utf    = boost::unit_test;
const double inf = std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(StrawSurfaces);
  /// Unit test for creating compliant/non-compliant StrawSurface object
  BOOST_AUTO_TEST_CASE(StrawSurfaceConstruction)
  {
    // StrawSurface default constructor is deleted
    //
    /// Constructor with transform pointer, null or valid, radius and halfZ
    double        radius(1.0), halfZ(10.);
    Translation3D translation{0., 1., 2.};
    auto pTransform     = std::make_shared<const Transform3D>(translation);
    auto pNullTransform = std::make_shared<const Transform3D>();
    BOOST_TEST(StrawSurface(pNullTransform, radius, halfZ).type()
               == Surface::Straw);
    BOOST_TEST(StrawSurface(pTransform, radius, halfZ).type()
               == Surface::Straw);
    //
    /// Constructor with transform and LineBounds pointer
    auto pLineBounds = std::make_shared<LineBounds>(radius, halfZ);
    BOOST_TEST(StrawSurface(pTransform, pLineBounds).type() == Surface::Straw);
    //
    /// Constructor with LineBounds ptr, DetectorElement and Identifier
    Identifier                                id{1};
    std::shared_ptr<const Acts::PlanarBounds> p
        = std::make_shared<const RectangleBounds>(1., 10.);
    DetectorElementStub detElement{id, pTransform, p, 1.0, nullptr};
    BOOST_TEST(StrawSurface(pLineBounds, detElement, id).type()
               == Surface::Straw);
    //
    /// Copy constructor
    StrawSurface strawSurfaceObject(pTransform, radius, halfZ);
    StrawSurface copiedStrawSurface(strawSurfaceObject);
    BOOST_TEST(copiedStrawSurface.type() == Surface::Straw);
    BOOST_TEST(copiedStrawSurface == strawSurfaceObject);
    //
    /// Copied and transformed
    StrawSurface copiedTransformedStrawSurface(strawSurfaceObject, *pTransform);
    BOOST_TEST(copiedTransformedStrawSurface.type() == Surface::Straw);
  }
  //
  /// Unit test for testing StrawSurface properties
  BOOST_AUTO_TEST_CASE(StrawSurfaceProperties)
  {
    /// Test clone method
    double        radius(1.0), halfZ(10.);
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    // auto pNullTransform = std::make_shared<const Transform3D>();
    StrawSurface strawSurfaceObject(pTransform, radius, halfZ);
    //
    auto pClonedStrawSurface = strawSurfaceObject.clone();
    BOOST_TEST(pClonedStrawSurface->type() == Surface::Straw);
    delete pClonedStrawSurface;
    //
    /// Test type (redundant)
    BOOST_TEST(strawSurfaceObject.type() == Surface::Straw);
    //
    /// Test name
    BOOST_TEST(strawSurfaceObject.name() == std::string("Acts::StrawSurface"));
    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    strawSurfaceObject.dump(dumpOuput);
    BOOST_TEST(dumpOuput.is_equal("Acts::StrawSurface\n\
     Center position  (x, y, z) = (0.0000, 1.0000, 2.0000)\n\
     Rotation:             colX = (1.000000, 0.000000, 0.000000)\n\
                           colY = (0.000000, 1.000000, 0.000000)\n\
                           colZ = (0.000000, 0.000000, 1.000000)\n\
     Bounds  : Acts::LineBounds: (radius, halflengthInZ) = (1.0000000, 10.0000000)"));
  }

  BOOST_AUTO_TEST_CASE(EqualityOperators, *utf::expected_failures(1))
  {
    double        radius(1.0), halfZ(10.);
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    StrawSurface  strawSurfaceObject(pTransform, radius, halfZ);
    //
    StrawSurface strawSurfaceObject2(pTransform, radius, halfZ);
    //
    /// Test equality operator
    BOOST_TEST(strawSurfaceObject == strawSurfaceObject2);
    //
    BOOST_TEST_CHECKPOINT(
        "Create and then assign a StrawSurface object to the existing one");
    /// Test assignment (will fail at the equality test)
    StrawSurface assignedStrawSurface(nullptr, NaN, NaN);
    assignedStrawSurface = strawSurfaceObject;
    /// Test equality of assigned to original
    BOOST_TEST(assignedStrawSurface == strawSurfaceObject);
  }
  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test

}  // end of namespace Acts
