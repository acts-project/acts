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
#include "Acts/Layers/Layer.hpp"
//#include "Acts/Utilities/Definitions.hpp"
#include "Acts/EventData/SingleTrackParameters.hpp"
#include "Acts/Layers/GenericApproachDescriptor.hpp"
#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include "Acts/Volumes/CuboidVolumeBounds.hpp"
#include "LayerStub.hpp"

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  namespace Layers {
    BOOST_AUTO_TEST_SUITE(Layers)

    /// Unit test for creating compliant/non-compliant Layer object
    BOOST_AUTO_TEST_CASE(LayerConstruction)
    {
      // Descendant Layer objects also inherit from Surface objects, which
      // delete the default constructor
      //
      /// Minimum possible construction (default constructor is deleted)
      LayerStub minallyConstructed(nullptr);
      BOOST_TEST(minallyConstructed.constructedOk());
      /// Need an approach descriptor for the next level of complexity:
      const std::vector<const Surface*> aSurfaces{new SurfaceStub(),
                                                  new SurfaceStub()};
      std::unique_ptr<ApproachDescriptor> ad(
          new GenericApproachDescriptor<Surface>(aSurfaces));
      const double thickness(1.0);
      LayerStub    approachDescriptorConstructed(
          nullptr, thickness, std::move(ad));
      /// Construction with (minimal) approach descriptor
      BOOST_TEST(approachDescriptorConstructed.constructedOk());
      // Copy construction is deleted
    }

    /// Unit test for testing Layer properties
    BOOST_AUTO_TEST_CASE(LayerProperties, *utf::expected_failures(1))
    {
      // Make a dummy layer to play with
      const std::vector<const Surface*> aSurfaces{new SurfaceStub(),
                                                  new SurfaceStub()};
      std::unique_ptr<ApproachDescriptor> ad(
          new GenericApproachDescriptor<Surface>(aSurfaces));
      auto                adPtr = ad.get();
      const double        thickness(1.0);
      SurfaceArrayCreator sac;
      double              halfX(0.1), halfY(0.2);
      size_t              binsX(2), binsY(4);
      auto                pSurfaceArray
          = sac.surfaceArrayOnPlane(aSurfaces, halfX, halfY, binsX, binsY);
      LayerStub layerStub(std::move(pSurfaceArray), thickness, std::move(ad));
      //
      /// surfaceArray()
      BOOST_TEST(layerStub.surfaceArray() == pSurfaceArray.get());
      /// thickness()
      BOOST_TEST(layerStub.thickness() == thickness);
      // onLayer() is templated; can't find implementation!
      /// isOnLayer() (delegates to the Surface 'isOnSurface()')
      const Vector3D pos{0.0, 0.0, 0.0};
      const Vector3D pos2{100., 100., std::nan("")};
      BOOST_TEST(layerStub.isOnLayer(pos) == true);
      // this should fail, but does not, but possibly my fault in SurfaceStub
      // implementation:
      BOOST_TEST(layerStub.isOnLayer(pos2) == false);
      /// approachDescriptor(), retrieved as a pointer.
      BOOST_TEST(layerStub.approachDescriptor() == adPtr);
      /// surfaceOnApproach()
      const Vector3D gpos{0., 0., 1.0};
      const Vector3D direction{0., 0., -1.};
      auto           surfaceOnApproach
          = layerStub.surfaceOnApproach(gpos,
                                        direction,
                                        PropDirection::alongMomentum,
                                        true,
                                        false,
                                        false,
                                        false);
      auto surfaceOnApproachIntersect = surfaceOnApproach.intersection;
      //(SurfaceStub uses hardcoded 20,20 dimensions)
      BOOST_TEST(surfaceOnApproachIntersect.pathLength == 20.0);
      // Need track parameters to test compatibleSurfaces
      SingleCurvilinearTrackParameters<ChargedPolicy> par(
          nullptr, {0., 0., 0.}, {0., 0., 1.}, -1);
      // Need more meaningful test, but docs need updating also; this is not the
      // signature in the doxygen documentation!
      auto compatibleSurfaces
          = layerStub.compatibleSurfaces(par,
                                         PropDirection::alongMomentum,
                                         true,
                                         false,
                                         false,
                                         false,
                                         0,
                                         nullptr,
                                         nullptr);
      BOOST_TEST(compatibleSurfaces.size() == 0u);
      /// nextLayer()
      BOOST_TEST(!(layerStub.nextLayer(gpos, direction)));
      /// enclosingTrackingVolume()
      BOOST_TEST(!layerStub.enclosingTrackingVolume());
      /// enclosingDetachedTrackingVolume()
      BOOST_TEST(!layerStub.enclosingDetachedTrackingVolume());
      /// registerRepresentingVolume(const AbstractVolume* avol)
      // need a volume:
      auto cubeVolumePtr = std::make_shared<CuboidVolumeBounds>(1., 2., 3.);
      AbstractVolume* abstractVolumePtr
          = new AbstractVolume(nullptr, cubeVolumePtr);
      layerStub.registerRepresentingVolume(abstractVolumePtr);
      BOOST_TEST(layerStub.representingVolume() == abstractVolumePtr);
      // BOOST_TEST_CHECKPOINT("Before ending test");
      // deletion results in "memory access violation at address: 0x00000071: no
      // mapping at fault address"
      // delete abstractVolumePtr;
      /// layerType()
      BOOST_TEST(layerStub.layerType() == LayerType::passive);
      /// detectorElements() (needs a better test)
      BOOST_TEST(layerStub.detectorElements().size() == 0u);
    }

    BOOST_AUTO_TEST_SUITE_END()
  }  // end of namespace Layers
}  // end of namespace Test

}  // end of namespace Acts
