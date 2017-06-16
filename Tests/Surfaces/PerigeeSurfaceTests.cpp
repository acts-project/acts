// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE PerigeeSurface Tests

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
#include "ACTS/Surfaces/InfiniteBounds.hpp"  //to get s_noBounds
#include "ACTS/Surfaces/PerigeeSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"  //to get s_noBounds
#include "ACTS/Utilities/Definitions.hpp"

using boost::test_tools::output_test_stream;
namespace utf    = boost::unit_test;
const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(PerigeeSurfaces);
  /// Unit test for creating compliant/non-compliant PerigeeSurface object
  BOOST_AUTO_TEST_CASE(PerigeeSurfaceConstruction)
  {
    // PerigeeSurface default constructor is deleted
    //
    /// Constructor with Vector3D
    Vector3D       unitXYZ{1., 1., 1.};
    PerigeeSurface perigeeSurfaceObject(unitXYZ);
    BOOST_TEST(PerigeeSurface(unitXYZ).type() == Surface::Perigee);
    //
    /// Constructor with transform pointer, null or valid
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    auto          pNullTransform = std::make_shared<const Transform3D>();
    BOOST_TEST(PerigeeSurface(pNullTransform).type() == Surface::Perigee);
    BOOST_TEST(PerigeeSurface(pTransform).type() == Surface::Perigee);
    //
    /// Copy constructor
    PerigeeSurface copiedPerigeeSurface(perigeeSurfaceObject);
    BOOST_TEST(copiedPerigeeSurface.type() == Surface::Perigee);
    BOOST_TEST(copiedPerigeeSurface == perigeeSurfaceObject);
    //
    /// Copied and transformed
    PerigeeSurface copiedTransformedPerigeeSurface(perigeeSurfaceObject,
                                                   *pTransform);
    BOOST_TEST(copiedTransformedPerigeeSurface.type() == Surface::Perigee);
  }
  //
  /// Unit test for testing PerigeeSurface properties
  BOOST_AUTO_TEST_CASE(PerigeeSurfaceProperties)
  {
    /// Test clone method
    Vector3D       unitXYZ{1., 1., 1.};
    PerigeeSurface perigeeSurfaceObject(unitXYZ);
    auto           pClonedPerigeeSurface = perigeeSurfaceObject.clone();
    BOOST_TEST(pClonedPerigeeSurface->type() == Surface::Perigee);
    delete pClonedPerigeeSurface;
    //
    /// Test type (redundant)
    BOOST_TEST(perigeeSurfaceObject.type() == Surface::Perigee);
    //
    /// Test name
    BOOST_TEST(perigeeSurfaceObject.name()
               == std::string("Acts::PerigeeSurface"));
    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    perigeeSurfaceObject.dump(dumpOuput);
    BOOST_TEST(dumpOuput.is_equal("Acts::PerigeeSurface:\n\
     Center position  (x, y, z) = (1.0000000, 1.0000000, 1.0000000)"));
  }

  BOOST_AUTO_TEST_CASE(EqualityOperators)
  {
    Vector3D       unitXYZ{1., 1., 1.};
    Vector3D       invalidPosition{NaN, NaN, NaN};
    PerigeeSurface perigeeSurfaceObject(unitXYZ);
    PerigeeSurface perigeeSurfaceObject2(unitXYZ);
    PerigeeSurface assignedPerigeeSurface(invalidPosition);
    /// Test equality operator
    BOOST_TEST(perigeeSurfaceObject == perigeeSurfaceObject2);
    /// Test assignment
    assignedPerigeeSurface = perigeeSurfaceObject;
    /// Test equality of assigned to original
    BOOST_TEST(assignedPerigeeSurface == perigeeSurfaceObject);
  }
  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test

}  // end of namespace Acts
