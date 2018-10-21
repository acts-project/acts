// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Layer Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/SingleTrackParameters.hpp"
#include "Acts/Layers/CylinderLayer.hpp"
#include "Acts/Layers/GenericApproachDescriptor.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include "Acts/Utilities/Definitions.hpp"
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
      BOOST_CHECK_EQUAL(pCylinderLayer->layerType(), LayerType::passive);
      // next level: need an array of Surfaces;
      // bounds object, rectangle type
      auto rBounds = std::make_shared<const RectangleBounds>(1., 1.);
      /// Constructor with transform pointer
      auto pNullTransform = std::make_shared<const Transform3D>();
      const std::vector<std::shared_ptr<const Surface>> aSurfaces{
          Surface::makeShared<PlaneSurface>(pNullTransform, rBounds),
          Surface::makeShared<PlaneSurface>(pNullTransform, rBounds)};
      const double thickness(1.0);
      auto         pCylinderLayerFromSurfaces
          = CylinderLayer::create(pTransform, pCylinder, nullptr);
      BOOST_CHECK_EQUAL(pCylinderLayerFromSurfaces->layerType(),
                        LayerType::passive);
      // construct with thickness:
      auto pCylinderLayerWithThickness
          = CylinderLayer::create(pTransform, pCylinder, nullptr, thickness);
      CHECK_CLOSE_REL(
          pCylinderLayerWithThickness->thickness(), thickness, 1e-6);
      // with an approach descriptor...
      std::unique_ptr<ApproachDescriptor> ad(
          new GenericApproachDescriptor(aSurfaces));
      auto adPtr                                = ad.get();
      auto pCylinderLayerWithApproachDescriptor = CylinderLayer::create(
          pTransform, pCylinder, nullptr, thickness, std::move(ad));
      BOOST_CHECK_EQUAL(
          pCylinderLayerWithApproachDescriptor->approachDescriptor(), adPtr);
      // with the layerType specified...
      auto pCylinderLayerWithLayerType
          = CylinderLayer::create(pTransform,
                                  pCylinder,
                                  nullptr,
                                  thickness,
                                  std::move(ad),
                                  LayerType::passive);
      BOOST_CHECK_EQUAL(pCylinderLayerWithLayerType->layerType(),
                        LayerType::passive);
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
      BOOST_CHECK_EQUAL(pCylinderLayer->surfaceRepresentation().name(),
                        std::string("Acts::CylinderSurface"));
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
      BOOST_CHECK_EQUAL(pl.get<double>("thickness"), 0.4);
      Transform3D act = from_variant<Transform3D>(pl.at("transform"));
      CHECK_CLOSE_OR_SMALL(*pTransform, act, 1e-6, 1e-9);

      auto pCylinderLayer2 = std::dynamic_pointer_cast<CylinderLayer>(
          CylinderLayer::create(var_data));

      BOOST_CHECK_EQUAL(pCylinderLayer->thickness(),
                        pCylinderLayer2->thickness());
      CHECK_CLOSE_OR_SMALL(pCylinderLayer->transform(),
                           pCylinderLayer2->transform(),
                           1e-6,
                           1e-9);

      auto cvBoundsExp = dynamic_cast<const CylinderVolumeBounds*>(
          &(pCylinderLayer->representingVolume()->volumeBounds()));
      auto cvBoundsAct = dynamic_cast<const CylinderVolumeBounds*>(
          &(pCylinderLayer2->representingVolume()->volumeBounds()));

      BOOST_CHECK_EQUAL(cvBoundsExp->innerRadius(), cvBoundsAct->innerRadius());
      BOOST_CHECK_EQUAL(cvBoundsExp->outerRadius(), cvBoundsAct->outerRadius());
      BOOST_CHECK_EQUAL(cvBoundsExp->halfPhiSector(),
                        cvBoundsAct->halfPhiSector());
      BOOST_CHECK_EQUAL(cvBoundsExp->halflengthZ(), cvBoundsAct->halflengthZ());
    }

    BOOST_AUTO_TEST_SUITE_END()
  }  // namespace Layers
}  // namespace Test

}  // namespace Acts
