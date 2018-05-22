// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
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
#include "Acts/Layers/NavigationLayer.hpp"
//#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include "Acts/Utilities/VariantData.hpp"
#include "Acts/Volumes/CuboidVolumeBounds.hpp"
#include "LayerStub.hpp"

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  namespace Layers {
    BOOST_AUTO_TEST_SUITE(Layers)

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

    BOOST_AUTO_TEST_CASE(NavigationLayer_toVariantData)
    {
      const double thickness = 0.1;
      double       w = 5, h = 10;
      auto         rbounds = std::make_shared<const RectangleBounds>(w, h);
      auto trf = std::make_shared<const Transform3D>(Translation3D(0, 0, 5));
      auto pSurface = std::make_unique<const PlaneSurface>(trf, rbounds);
      auto pNavigationLayer = std::dynamic_pointer_cast<const NavigationLayer>(
          NavigationLayer::create(std::move(pSurface), thickness));

      variant_data var_data = pNavigationLayer->toVariantData();
      std::cout << var_data << std::endl;

      auto pNavigationLayer2 = std::dynamic_pointer_cast<const NavigationLayer>(
          NavigationLayer::create(var_data));
      std::cout << pNavigationLayer2->toVariantData() << std::endl;

      auto rbounds2 = dynamic_cast<const RectangleBounds*>(
          &pNavigationLayer2->surfaceRepresentation().bounds());
      BOOST_TEST(trf->isApprox(
          pNavigationLayer2->surfaceRepresentation().transform()));
      BOOST_TEST(rbounds2->halflengthX() == w);
      BOOST_TEST(rbounds2->halflengthY() == h);
    }

    BOOST_AUTO_TEST_SUITE_END()
  }  // end of namespace Layers
}  // end of namespace Test

}  // end of namespace Acts
