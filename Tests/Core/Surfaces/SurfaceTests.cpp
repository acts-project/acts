// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Surface Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Layers/PlaneLayer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"   //to get s_noBounds
#include "Acts/Surfaces/RectangleBounds.hpp"  //to get s_noBounds
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include "SurfaceStub.hpp"

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {
/// Mock track object with minimal methods implemented for compilation
class MockTrack
{
public:
  MockTrack(const Vector3D& mom, const Vector3D& pos) : m_mom(mom), m_pos(pos)
  {
    // nop
  }

  Vector3D
  momentum() const
  {
    return m_mom;
  }

  Vector3D
  position() const
  {
    return m_pos;
  }

private:
  Vector3D m_mom;
  Vector3D m_pos;
};

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)

  /// todo: make test fixture; separate out different cases

  /// Unit test for creating compliant/non-compliant Surface object
  BOOST_AUTO_TEST_CASE(SurfaceConstruction)
  {
    // SurfaceStub s;
    BOOST_CHECK_EQUAL(Surface::Other, SurfaceStub().type());
    SurfaceStub original;
    BOOST_CHECK_EQUAL(Surface::Other, SurfaceStub(original).type());
    Translation3D translation{0., 1., 2.};
    Transform3D   transform(translation);
    BOOST_CHECK_EQUAL(Surface::Other, SurfaceStub(original, transform).type());
    // need some cruft to make the next one work
    auto pTransform = std::make_shared<const Transform3D>(translation);
    std::shared_ptr<const Acts::PlanarBounds> p
        = std::make_shared<const RectangleBounds>(5., 10.);
    DetectorElementStub detElement{pTransform, p, 0.2, nullptr};
    BOOST_CHECK_EQUAL(Surface::Other, SurfaceStub(detElement).type());
  }

  /// Unit test for testing Surface properties
  BOOST_AUTO_TEST_CASE(SurfaceProperties, *utf::expected_failures(1))
  {
    // build a test object , 'surface'
    std::shared_ptr<const Acts::PlanarBounds> pPlanarBound
        = std::make_shared<const RectangleBounds>(5., 10.);
    Vector3D      reference{0., 1., 2.};
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    auto          pLayer     = PlaneLayer::create(pTransform, pPlanarBound);
    MaterialProperties properties{0.2, 0.2, 0.2, 20., 10, 5.};
    auto               pMaterial
        = std::make_shared<const HomogeneousSurfaceMaterial>(properties);
    DetectorElementStub detElement{pTransform, pPlanarBound, 0.2, pMaterial};
    SurfaceStub         surface(detElement);
    // associatedDetectorElement
    BOOST_CHECK_EQUAL(surface.associatedDetectorElement(), &detElement);
    // test  associatelayer, associatedLayer
    surface.associateLayer(*pLayer);
    BOOST_CHECK_EQUAL(surface.associatedLayer(), pLayer.get());
    // associated Material is not set to the surface
    // it is set to the detector element surface though
    BOOST_CHECK_NE(surface.associatedMaterial(), pMaterial.get());
    // center()
    CHECK_CLOSE_OR_SMALL(reference, surface.center(), 1e-6, 1e-9);
    // stream insertion operator <<
    output_test_stream output;
    output << surface;
    BOOST_CHECK(!output.is_empty(false));  // no check on contents
    // insideBounds
    Vector2D localPosition{0.1, 3.0};
    BOOST_CHECK(surface.insideBounds(localPosition));
    Vector2D outside{20., 20.};
    BOOST_CHECK(!surface.insideBounds(
        outside));  // fails: m_bounds only in derived classes
    // intersectionEstimate (should delegate to derived class method of same
    // name)
    Vector3D mom{100., 200., 300.};
    auto     intersectionEstimate
        = surface.intersectionEstimate(reference, mom, forward, false);
    const Intersection ref{Vector3D{1, 1, 1}, 20., true};
    BOOST_CHECK_EQUAL(ref.position, intersectionEstimate.position);
    // isFree
    BOOST_CHECK(!surface.isFree());
    // isOnSurface
    BOOST_CHECK(surface.isOnSurface(reference, mom, false));
    BOOST_CHECK(
        surface.isOnSurface(reference, mom, true));  // need to improve bounds()
    // referenceFrame()
    RotationMatrix3D unitary;
    unitary << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    auto referenceFrame = surface.referenceFrame(
        reference, mom);  // need more complex case to test
    BOOST_CHECK_EQUAL(referenceFrame, unitary);
    // normal()
    auto normal = surface.Surface::normal(reference);  // needs more complex
                                                       // test
    Vector3D zero{0., 0., 0.};
    BOOST_CHECK_EQUAL(zero, normal);
    // pathCorrection is pure virtual
    // associatedMaterial()
    MaterialProperties newProperties{0.5, 0.5, 0.5, 20., 10., 5.};
    auto               pNewMaterial
        = std::make_shared<const HomogeneousSurfaceMaterial>(newProperties);
    surface.setAssociatedMaterial(pNewMaterial);
    BOOST_CHECK_EQUAL(surface.associatedMaterial(),
                      pNewMaterial.get());  // passes ??
    //
    CHECK_CLOSE_OR_SMALL(surface.transform(), *pTransform, 1e-6, 1e-9);
    // type() is pure virtual
  }

  BOOST_AUTO_TEST_CASE(EqualityOperators)
  {
    // build some test objects
    std::shared_ptr<const Acts::PlanarBounds> pPlanarBound
        = std::make_shared<const RectangleBounds>(5., 10.);
    Vector3D      reference{0., 1., 2.};
    Translation3D translation1{0., 1., 2.};
    Translation3D translation2{1., 1., 2.};
    auto pTransform1 = std::make_shared<const Transform3D>(translation1);
    auto pTransform2 = std::make_shared<const Transform3D>(translation2);
    auto pLayer      = PlaneLayer::create(pTransform1, pPlanarBound);
    MaterialProperties properties{1., 1., 1., 20., 10, 5.};
    auto               pMaterial
        = std::make_shared<const HomogeneousSurfaceMaterial>(properties);
    DetectorElementStub detElement1{pTransform1, pPlanarBound, 0.2, pMaterial};
    DetectorElementStub detElement2{pTransform1, pPlanarBound, 0.3, pMaterial};
    DetectorElementStub detElement3{pTransform2, pPlanarBound, 0.3, pMaterial};
    //
    SurfaceStub surface1(detElement1);
    SurfaceStub surface2(detElement1);  // 1 and 2 are the same
    SurfaceStub surface3(detElement2);  // 3 differs in thickness
    SurfaceStub surface4(detElement3);  // 4 has a different transform and id
    //
    BOOST_CHECK_EQUAL(surface1, surface2);
    //
    // remove test for the moment,
    // surfaces do not have a concept of thickness (only detector elements have)
    // only thickness is different here
    //
    // BOOST_CHECK_NE(surface1, surface3);  // will fail
    //
    BOOST_CHECK_NE(surface1, surface4);
  }
  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
