// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Layer Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//#include <limits>
#include "ACTS/Layers/NavigationLayer.hpp"
//#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Volumes/CuboidVolumeBounds.hpp"
#include "LayerStub.hpp"

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  namespace Layers {
    BOOST_AUTO_TEST_SUITE(Layers);

    /// Unit test for creating compliant/non-compliant NavigationLayer object
    BOOST_AUTO_TEST_CASE(NavigationLayerConstruction)
    {
      // default constructor, copy and assignment are all deleted
      auto pSurface         = std::unique_ptr<const Surface>(new SurfaceStub());
      auto pNavigationLayer = NavigationLayer::create(std::move(pSurface));
      BOOST_TEST(pNavigationLayer->layerType() == LayerType::navigation);
      // next level: with thickness
      const double thickness = 0.1;
      auto pSurface2 = std::unique_ptr<const Surface>(new SurfaceStub());
      auto pThickNavigationLayer
          = NavigationLayer::create(std::move(pSurface2), thickness);
      BOOST_TEST(pThickNavigationLayer->layerType() == LayerType::navigation);
    }

    /// Unit test for testing NavigationLayer properties
    BOOST_AUTO_TEST_CASE(NavigationLayerProperties, *utf::expected_failures(1))
    {
      const double thickness = 0.1;
      auto         pSurface = std::unique_ptr<const Surface>(new SurfaceStub());
      auto         rawSurfacePtr = pSurface.get();
      auto         pNavigationLayer
          = NavigationLayer::create(std::move(pSurface), thickness);
      BinningValue b{BinningValue::binZ};
      Vector3D     origin{0., 0., 0.};
      // binningPosition(), needs a better test
      BOOST_TEST(pNavigationLayer->binningPosition(b) == origin);
      // surfaceRepresentation() [looks dangerous]
      BOOST_TEST(rawSurfacePtr == &(pNavigationLayer->surfaceRepresentation()));
      // isOnLayer()
      BOOST_CHECK(pNavigationLayer->isOnLayer(origin, true));
      // isOnLayer()
      Vector3D crazyPosition{1000., 10000., std::nan("")};
      BOOST_CHECK(pNavigationLayer->isOnLayer(crazyPosition, true) == false);
      // resolve()
      BOOST_CHECK(pNavigationLayer->resolve(true, true, true) == false);
    }

    BOOST_AUTO_TEST_SUITE_END();
  }  // end of namespace Layers
}  // end of namespace Test

}  // end of namespace Acts
