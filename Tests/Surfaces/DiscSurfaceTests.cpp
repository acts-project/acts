// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Disc Surface Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//
#include "ACTS/Material/HomogeneousSurfaceMaterial.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
//
#include "DetectorElementStub.hpp"
//
#include <limits>

namespace utf    = boost::unit_test;
const double inf = std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  // using boost::test_tools::output_test_stream;
  
  BOOST_AUTO_TEST_SUITE(Surfaces);
  /// Unit tests for creating DiscSurface object
  BOOST_AUTO_TEST_CASE(DiscSurface_constructors_test)
  {
    //default constructor is deleted
    //scaffolding...
    double rMin(1.0), rMax(5.0), halfPhiSector(M_PI/8.);
    //
    ///Test DiscSurface fully specified constructor but no transform
    BOOST_CHECK_NO_THROW(DiscSurface(nullptr, rMin, rMax, halfPhiSector));
    //
    ///Test DiscSurface constructor with default halfPhiSector
    BOOST_CHECK_NO_THROW(DiscSurface(nullptr, rMin, rMax));
    //
    ///Test DiscSurface constructor with a transform specified
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<Transform3D>(translation);
    BOOST_CHECK_NO_THROW(DiscSurface(pTransform, rMin, rMax, halfPhiSector));
    //
    ///Copy constructed DiscSurface
    DiscSurface anotherDiscSurface(pTransform, rMin, rMax, halfPhiSector);
    //N.B. Just using BOOST_CHECK_NO_THROW(DiscSurface(anotherDiscSurface)) tries to call 
    //the (deleted) default constructor.
    DiscSurface copiedDiscSurface(anotherDiscSurface);
    BOOST_TEST_MESSAGE("Copy constructed DiscSurface ok");
    //
    ///Copied and transformed DiscSurface
    BOOST_CHECK_NO_THROW(DiscSurface(anotherDiscSurface, *pTransform));
  }
  
  /// Unit tests of all named methods
  BOOST_AUTO_TEST_CASE(DiscSurface_properties_test)
  {
    Vector3D origin3D{0,0,0};
    std::shared_ptr<Transform3D> pTransform; //nullptr
    double rMin(1.0), rMax(1.0), halfPhiSector(M_PI/8.);
    DiscSurface discSurfaceObject(pTransform, rMin, rMax, halfPhiSector);
    //
    ///Test type
    BOOST_CHECK(discSurfaceObject.type() == Surface::Disc);
    //
    ///Test normal, no local position specified
    Vector3D zAxis{0,0,1};
    BOOST_CHECK(discSurfaceObject.normal() == zAxis);
    //
    ///Test normal, local position specified
    Vector2D lpos(2.0,0.05);
    BOOST_CHECK(discSurfaceObject.normal(lpos) == zAxis);
    //
    ///Test binningPosition
    auto binningPosition= discSurfaceObject.binningPosition(BinningValue::binRPhi );
    std::cout<<binningPosition<<std::endl;
    //
    ///Test bounds 
    BOOST_TEST(discSurfaceObject.bounds().type() = SurfaceBounds::Disc);
    //
    ///Test isOnSurface()
    BOOST_TEST(discSurfaceObject.isOnSurface(origin3D, true));
  }
  /// Unit test for testing DiscSurface assignment
  BOOST_AUTO_TEST_CASE(DiscSurface_assignment_test)
  {
    
  }
  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test

}  // end of namespace Acts
