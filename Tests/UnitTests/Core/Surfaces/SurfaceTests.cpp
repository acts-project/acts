// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include <limits>

#include "Acts/Geometry/PlaneLayer.hpp"
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
class MockTrack {
 public:
  MockTrack(const Vector3D& mom, const Vector3D& pos) : m_mom(mom), m_pos(pos) {
    // nop
  }

  Vector3D momentum() const { return m_mom; }

  Vector3D position() const { return m_pos; }

 private:
  Vector3D m_mom;
  Vector3D m_pos;
};

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(Surfaces)

/// todo: make test fixture; separate out different cases

/// Unit test for creating compliant/non-compliant Surface object
BOOST_AUTO_TEST_CASE(SurfaceConstruction) {
  // SurfaceStub s;
  BOOST_CHECK_EQUAL(Surface::Other, SurfaceStub().type());
  SurfaceStub original;
  BOOST_CHECK_EQUAL(Surface::Other, SurfaceStub(original).type());
  Translation3D translation{0., 1., 2.};
  Transform3D transform(translation);
  BOOST_CHECK_EQUAL(Surface::Other,
                    SurfaceStub(tgContext, original, transform).type());
  // need some cruft to make the next one work
  auto pTransform = std::make_shared<const Transform3D>(translation);
  std::shared_ptr<const Acts::PlanarBounds> p =
      std::make_shared<const RectangleBounds>(5., 10.);
  DetectorElementStub detElement{pTransform, p, 0.2, nullptr};
  BOOST_CHECK_EQUAL(Surface::Other, SurfaceStub(detElement).type());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
