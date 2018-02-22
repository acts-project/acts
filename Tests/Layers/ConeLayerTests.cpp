// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
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
#include "ACTS/Layers/ConeLayer.hpp"
#include "ACTS/Surfaces/ConeBounds.hpp"
//#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/EventData/SingleTrackParameters.hpp"
#include "ACTS/Layers/GenericApproachDescriptor.hpp"
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Volumes/CuboidVolumeBounds.hpp"
#include "LayerStub.hpp"

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  namespace Layers {
    BOOST_AUTO_TEST_SUITE(Layers)

    /// Unit test for creating compliant/non-compliant ConeLayer object
    BOOST_AUTO_TEST_CASE(ConeLayerConstruction)
    {
      // default constructor, copy and assignment are all deleted
      // minimally need a Transform3D and a PlanarBounds object (e.g.
      // ConeBounds) to construct
      Translation3D translation{0., 1., 2.};
      auto       pTransform = std::make_shared<const Transform3D>(translation);
      double     alpha(M_PI / 8.0);
      const bool symmetric(false);
      auto       pCone = std::make_shared<const ConeBounds>(alpha, symmetric);
      // for some reason, this one doesnt exist
      // auto         pConeLayer = ConeLayer::create(pTransform, pCone);
      // BOOST_TEST(pConeLayer->layerType() == LayerType::passive);
      // next level: need an array of Surfaces;
      const std::vector<const Surface*> aSurfaces{new SurfaceStub(),
                                                  new SurfaceStub()};
      const double        thickness(1.0);
      SurfaceArrayCreator sac;
      size_t              binsX(2), binsY(4);
      auto                pSurfaceArray
          = sac.surfaceArrayOnPlane(aSurfaces, 10, 20, binsX, binsY);
      auto pConeLayerFromSurfaces
          = ConeLayer::create(pTransform, pCone, std::move(pSurfaceArray));
      BOOST_TEST(pConeLayerFromSurfaces->layerType() == LayerType::active);
      // construct with thickness:
      auto pConeLayerWithThickness = ConeLayer::create(
          pTransform, pCone, std::move(pSurfaceArray), thickness);
      BOOST_TEST(pConeLayerWithThickness->thickness() == thickness);
      // with an approach descriptor...
      std::unique_ptr<ApproachDescriptor> ad(
          new GenericApproachDescriptor<Surface>(aSurfaces));
      auto adPtr = ad.get();
      auto pConeLayerWithApproachDescriptor
          = ConeLayer::create(pTransform,
                              pCone,
                              std::move(pSurfaceArray),
                              thickness,
                              std::move(ad));
      BOOST_TEST(pConeLayerWithApproachDescriptor->approachDescriptor()
                 == adPtr);
      // with the layerType specified...
      auto pConeLayerWithLayerType = ConeLayer::create(pTransform,
                                                       pCone,
                                                       std::move(pSurfaceArray),
                                                       thickness,
                                                       std::move(ad),
                                                       LayerType::passive);
      BOOST_TEST(pConeLayerWithLayerType->layerType() == LayerType::passive);
    }

    /// Unit test for testing Layer properties
    BOOST_AUTO_TEST_CASE(ConeLayerProperties /*, *utf::expected_failures(1)*/)
    {
      Translation3D translation{0., 1., 2.};
      auto       pTransform = std::make_shared<const Transform3D>(translation);
      double     alpha(M_PI / 8.0);
      const bool symmetric(false);
      auto       pCone = std::make_shared<const ConeBounds>(alpha, symmetric);
      const std::vector<const Surface*> aSurfaces{new SurfaceStub(),
                                                  new SurfaceStub()};
      // const double        thickness(1.0);
      SurfaceArrayCreator sac;
      size_t              binsX(2), binsY(4);
      auto                pSurfaceArray
          = sac.surfaceArrayOnPlane(aSurfaces, 10, 20, binsX, binsY);
      auto pConeLayerFromSurfaces
          = ConeLayer::create(pTransform, pCone, std::move(pSurfaceArray));
      // auto planeSurface = pConeLayer->surfaceRepresentation();
      BOOST_TEST(pConeLayerFromSurfaces->surfaceRepresentation().name()
                 == std::string("Acts::ConeSurface"));
    }

    BOOST_AUTO_TEST_SUITE_END()
  }  // end of namespace Layers
}  // end of namespace Test

}  // end of namespace Acts
