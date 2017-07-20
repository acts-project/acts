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

#include "../Surfaces/SurfaceStub.hpp"
#include "ACTS/Layers/GenericApproachDescriptor.hpp"
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "LayerStub.hpp"

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  namespace Layers {
    BOOST_AUTO_TEST_SUITE(Layers);

    /// Unit test for creating compliant/non-compliant GenericApproachDescriptor
    /// object
    BOOST_AUTO_TEST_CASE(GenericApproachDescriptorConstruction)
    {
      typedef std::vector<const Surface*> SurfaceVector_t;
      //
      const SurfaceVector_t someSurfaces{new SurfaceStub(), new SurfaceStub()};
      BOOST_CHECK_NO_THROW(
          GenericApproachDescriptor<Surface>
              minimallyConstructedApproachDescriptor(someSurfaces));
      //
      std::vector<std::shared_ptr<const Layer>> sharedLayers{
          std::make_shared<LayerStub>(nullptr),
          std::make_shared<LayerStub>(nullptr)};
      BOOST_CHECK_NO_THROW(
          GenericApproachDescriptor<Layer> sharedLayerApproachDescriptor(
              sharedLayers));
    }

    /// Unit test for testing GenericApproachDescriptor properties
    BOOST_AUTO_TEST_CASE(GenericApproachDescriptorProperties,
                         *utf::expected_failures(1))
    {
      typedef std::vector<const Surface*> SurfaceVector_t;
      Vector3D                            origin{
          0., 0., 0.,
      };
      Vector3D      zDir{0., 0., 1.};
      BoundaryCheck bcheck{true};
      //
      const SurfaceVector_t someSurfaces{new SurfaceStub(), new SurfaceStub()};
      GenericApproachDescriptor<Surface> approachDescriptor(someSurfaces);
      LayerStub                          aLayer(nullptr);
      // registerLayer()
      BOOST_CHECK_NO_THROW(approachDescriptor.registerLayer(aLayer));
      // approachSurface
      SurfaceIntersection surfIntersection
          = approachDescriptor.approachSurface(origin, zDir, bcheck);
      double expectedIntersection = 20.0;  // property of SurfaceStub
      BOOST_CHECK(surfIntersection.intersection.pathLength
                  == expectedIntersection);
      // containedSurfaces()
      BOOST_CHECK(&(approachDescriptor.containedSurfaces()) == &someSurfaces);
    }

    BOOST_AUTO_TEST_SUITE_END();
  }  // end of namespace Layers
}  // end of namespace Test

}  // end of namespace Acts
