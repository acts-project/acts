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
#include "Acts/Layers/CylinderLayer.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
//#include "Acts/Utilities/Definitions.hpp"
#include "Acts/EventData/SingleTrackParameters.hpp"
#include "Acts/Layers/GenericApproachDescriptor.hpp"
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

    /// Unit test for creating compliant/non-compliant CylinderLayer object
    BOOST_AUTO_TEST_CASE(CylinderLayerConstruction)
    {
      // default constructor, copy and assignment are all deleted
      // minimally need a Transform3D and a PlanarBounds object (e.g.
      // CylinderBounds) to construct
      Translation3D translation{0., 1., 2.};
      auto   pTransform = std::make_shared<const Transform3D>(translation);
      double radius(0.5), halfz(10.);
      auto   pCylinder = std::make_shared<const CylinderBounds>(radius, halfz);
      auto   pCylinderLayer = CylinderLayer::create(pTransform, pCylinder);
      BOOST_TEST(pCylinderLayer->layerType() == LayerType::passive);
      // next level: need an array of Surfaces;
      const std::vector<const Surface*> aSurfaces{new SurfaceStub(),
                                                  new SurfaceStub()};
      const double        thickness(1.0);
      SurfaceArrayCreator sac;
      size_t              binsX(2), binsY(4);
      auto                pSurfaceArray
          = sac.surfaceArrayOnPlane(aSurfaces, 10, 20, binsX, binsY);
      auto pCylinderLayerFromSurfaces = CylinderLayer::create(
          pTransform, pCylinder, std::move(pSurfaceArray));
      BOOST_TEST(pCylinderLayerFromSurfaces->layerType() == LayerType::passive);
      // construct with thickness:
      auto pCylinderLayerWithThickness = CylinderLayer::create(
          pTransform, pCylinder, std::move(pSurfaceArray), thickness);
      BOOST_CHECK_CLOSE_FRACTION(
          pCylinderLayerWithThickness->thickness(), thickness, 1e-6);
      // with an approach descriptor...
      std::unique_ptr<ApproachDescriptor> ad(
          new GenericApproachDescriptor<Surface>(aSurfaces));
      auto adPtr = ad.get();
      auto pCylinderLayerWithApproachDescriptor
          = CylinderLayer::create(pTransform,
                                  pCylinder,
                                  std::move(pSurfaceArray),
                                  thickness,
                                  std::move(ad));
      BOOST_TEST(pCylinderLayerWithApproachDescriptor->approachDescriptor()
                 == adPtr);
      // with the layerType specified...
      auto pCylinderLayerWithLayerType
          = CylinderLayer::create(pTransform,
                                  pCylinder,
                                  std::move(pSurfaceArray),
                                  thickness,
                                  std::move(ad),
                                  LayerType::passive);
      BOOST_TEST(pCylinderLayerWithLayerType->layerType()
                 == LayerType::passive);
    }

    /// Unit test for testing Layer properties
    BOOST_AUTO_TEST_CASE(
        CylinderLayerProperties /*, *utf::expected_failures(1)*/)
    {
      Translation3D translation{0., 1., 2.};
      auto   pTransform = std::make_shared<const Transform3D>(translation);
      double radius(0.5), halfz(10.);
      auto   pCylinder = std::make_shared<const CylinderBounds>(radius, halfz);
      auto   pCylinderLayer = CylinderLayer::create(pTransform, pCylinder);
      // auto planeSurface = pCylinderLayer->surfaceRepresentation();
      BOOST_TEST(pCylinderLayer->surfaceRepresentation().name()
                 == std::string("Acts::CylinderSurface"));
    }

    BOOST_AUTO_TEST_CASE(CylinderLayer_toVariantData)
    {
      Translation3D translation{0., 1., 2.};
      Transform3D   rot;
      rot = AngleAxis3D(M_PI / 4, Vector3D::UnitZ());

      auto pTransform = std::make_shared<const Transform3D>(translation * rot);
      double radius(0.5), halfz(10.);
      auto   pCylinder = std::make_shared<const CylinderBounds>(radius, halfz);
      auto   pCylinderLayer = std::dynamic_pointer_cast<CylinderLayer>(
          CylinderLayer::create(pTransform, pCylinder, nullptr, 0.4));

      variant_data var_data = pCylinderLayer->toVariantData();
      std::cout << var_data << std::endl;

      variant_map var_map = boost::get<variant_map>(var_data);
      variant_map pl      = var_map.get<variant_map>("payload");
      BOOST_CHECK_CLOSE_FRACTION(pl.get<double>("thickness"), 0.4, 1e-6);
      Transform3D act = from_variant<Transform3D>(pl.at("transform"));
      BOOST_TEST(pTransform->isApprox(act));

      auto pCylinderLayer2 = std::dynamic_pointer_cast<CylinderLayer>(
          CylinderLayer::create(var_data));

      BOOST_CHECK_CLOSE_FRACTION(
          pCylinderLayer->thickness(), pCylinderLayer2->thickness(), 1e-6);
      BOOST_TEST(
          pCylinderLayer->transform().isApprox(pCylinderLayer2->transform()));

      auto cvBoundsExp = dynamic_cast<const CylinderVolumeBounds*>(
          &(pCylinderLayer->representingVolume()->volumeBounds()));
      auto cvBoundsAct = dynamic_cast<const CylinderVolumeBounds*>(
          &(pCylinderLayer2->representingVolume()->volumeBounds()));

      BOOST_CHECK_CLOSE_FRACTION(
          cvBoundsExp->innerRadius(), cvBoundsAct->innerRadius(), 1e-6);
      BOOST_CHECK_CLOSE_FRACTION(
          cvBoundsExp->outerRadius(), cvBoundsAct->outerRadius(), 1e-6);
      BOOST_CHECK_CLOSE_FRACTION(
          cvBoundsExp->halfPhiSector(), cvBoundsAct->halfPhiSector(), 1e-6);
      BOOST_CHECK_CLOSE_FRACTION(
          cvBoundsExp->halflengthZ(), cvBoundsAct->halflengthZ(), 1e-6);
    }

    BOOST_AUTO_TEST_SUITE_END()
  }  // namespace Layers
}  // namespace Test

}  // namespace Acts
