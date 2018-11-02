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

#include "../Surfaces/SurfaceStub.hpp"
#include "Acts/Layers/GenericApproachDescriptor.hpp"
#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include "LayerStub.hpp"

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  namespace Layers {
    BOOST_AUTO_TEST_SUITE(Layers)

    /// Unit test for creating compliant/non-compliant GenericApproachDescriptor
    /// object
    BOOST_AUTO_TEST_CASE(GenericApproachDescriptorConstruction)
    {
      std::vector<std::shared_ptr<const Surface>> someSurfaces{
          Surface::makeShared<SurfaceStub>(),
          Surface::makeShared<SurfaceStub>()};
      BOOST_CHECK_NO_THROW(
          GenericApproachDescriptor minimallyConstructedApproachDescriptor(
              someSurfaces));
      //
      std::vector<std::shared_ptr<const Layer>> sharedLayers{
          std::make_shared<LayerStub>(nullptr),
          std::make_shared<LayerStub>(nullptr)};
      BOOST_CHECK_NO_THROW(
          GenericApproachDescriptor sharedLayerApproachDescriptor(
              {sharedLayers.at(0)->surfaceRepresentation().getSharedPtr(),
               sharedLayers.at(1)->surfaceRepresentation().getSharedPtr()}));
    }

    /// Unit test for testing GenericApproachDescriptor properties
    BOOST_AUTO_TEST_CASE(GenericApproachDescriptorProperties,
                         *utf::expected_failures(1))
    {
      Vector3D origin{
          0., 0., 0.,
      };
      Vector3D      zDir{0., 0., 1.};
      BoundaryCheck bcheck{true};
      //
      std::vector<std::shared_ptr<const Surface>> someSurfaces{
          Surface::makeShared<SurfaceStub>(),
          Surface::makeShared<SurfaceStub>()};
      GenericApproachDescriptor approachDescriptor(someSurfaces);
      LayerStub                 aLayer(nullptr);
      // registerLayer()
      BOOST_CHECK_NO_THROW(approachDescriptor.registerLayer(aLayer));
      // approachSurface
      SurfaceIntersection surfIntersection
          = approachDescriptor.approachSurface(origin, zDir, forward, bcheck);
      double expectedIntersection = 20.0;  // property of SurfaceStub
      BOOST_CHECK(surfIntersection.intersection.pathLength
                  == expectedIntersection);
      // containedSurfaces()
      BOOST_TEST(approachDescriptor.containedSurfaces().size()
                 == someSurfaces.size());

      for (size_t i = 0; i < someSurfaces.size(); i++) {
        BOOST_TEST(approachDescriptor.containedSurfaces().at(i)
                   == someSurfaces.at(i).get());
      }
      // BOOST_CHECK(&(approachDescriptor.containedSurfaces()) ==
      // &someSurfaces);
    }

    BOOST_AUTO_TEST_SUITE_END()
  }  // namespace Layers
}  // namespace Test

}  // namespace Acts
