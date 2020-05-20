// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/ConeVolumeBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/ObjTestWriter.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace tt = boost::test_tools;

namespace Acts {

using namespace UnitLiterals;

// Create a test context
GeometryContext tgContext = GeometryContext();

namespace Test {

BOOST_AUTO_TEST_SUITE(VolumeBounds)

BOOST_AUTO_TEST_CASE(ConeVolumeBoundsTests) {
  std::vector<IdentifiedPolyhedron> tPolyhedrons;

  // Single solid Cone
  ConeVolumeBounds solidCone(0., 0., 0.45, 50_mm, 50_mm, 0., M_PI);

  // Test correct parameter return
  BOOST_TEST(solidCone.get(ConeVolumeBounds::eInnerAlpha) == 0.);
  BOOST_TEST(solidCone.get(ConeVolumeBounds::eInnerOffsetZ) == 0.);
  BOOST_TEST(solidCone.get(ConeVolumeBounds::eOuterAlpha) == 0.45);
  BOOST_TEST(solidCone.get(ConeVolumeBounds::eOuterOffsetZ) == 50.);
  BOOST_TEST(solidCone.get(ConeVolumeBounds::eHalfLengthZ) == 50.);
  BOOST_TEST(solidCone.get(ConeVolumeBounds::eAveragePhi) == 0.);
  BOOST_TEST(solidCone.get(ConeVolumeBounds::eHalfPhiSector) == M_PI);
  // Derived quantities
  BOOST_TEST(solidCone.innerTanAlpha() == 0.);
  BOOST_TEST(solidCone.innerRmin() == 0.);
  BOOST_TEST(solidCone.innerRmax() == 0.);
  BOOST_TEST(solidCone.outerTanAlpha() == std::tan(0.45));

  double outerRmax = 100_mm * solidCone.outerTanAlpha();
  BOOST_TEST(solidCone.outerRmin() == 0.);
  BOOST_TEST(solidCone.outerRmax() == outerRmax);

  auto solidConeSurfaces = solidCone.orientedSurfaces();
  BOOST_TEST(solidConeSurfaces.size() == 2);

  // Single solid Cone - with cut off
  ConeVolumeBounds cutOffCone(0., 0., 0.45, 80_mm, 50_mm, 0., M_PI);
  auto cutOffConeSurfaces = cutOffCone.orientedSurfaces();
  BOOST_TEST(cutOffConeSurfaces.size() == 3);

  // Cone - Cone inlay
  ConeVolumeBounds cutOffHollowCone(0.35, 70_mm, 0.45, 80_mm, 50_mm, 0., M_PI);
  auto cutOffHollowConeSurfaces = cutOffHollowCone.orientedSurfaces();
  BOOST_TEST(cutOffHollowConeSurfaces.size() == 4);

  // Sectoral Cone - Cone inlay
  ConeVolumeBounds cutOffHollowSectoralCone(0.35, 70_mm, 0.45, 80_mm, 50_mm, 0.,
                                            0.456);
  auto cutOffHollowSectoralConeSurfaces =
      cutOffHollowSectoralCone.orientedSurfaces();
  BOOST_TEST(cutOffHollowSectoralConeSurfaces.size() == 6);

  // Sectoral Cone - Hollow Cone
  ConeVolumeBounds cutOffHollowCylCone(10_mm, 0.45, 80_mm, 50_mm, 0., M_PI);
  auto cutOffHollowCylConeSurfaces = cutOffHollowCylCone.orientedSurfaces();
  BOOST_TEST(cutOffHollowCylConeSurfaces.size() == 4);

  // Single Hollow Cylinder - Cone inlay
  ConeVolumeBounds cutOffHollowConeCyl(120_mm, 0.35, 70_mm, 50_mm, 0., M_PI);
  auto cutOffHollowConeCylSurfaces = cutOffHollowConeCyl.orientedSurfaces();
  BOOST_TEST(cutOffHollowConeCylSurfaces.size() == 4);
}

BOOST_AUTO_TEST_CASE(ConeVolumeBoundsSurfaceOrientation) {
  GeometryContext tgContext = GeometryContext();

  ConeVolumeBounds hcone(10_mm, 0.45, 80_mm, 50_mm, 0., M_PI);

  auto cvbOrientedSurfaces = hcone.orientedSurfaces(nullptr);
  BOOST_TEST(cvbOrientedSurfaces.size(), 4);

  auto geoCtx = GeometryContext();
  Vector3D xaxis(1., 0., 0.);
  Vector3D yaxis(0., 1., 0.);
  Vector3D zaxis(0., 0., 1.);

  for (auto& os : cvbOrientedSurfaces) {
    // Test the orientation of the boundary surfaces
    auto rot = os.first->transform(geoCtx).rotation();
    BOOST_CHECK(rot.col(0).isApprox(xaxis));
    BOOST_CHECK(rot.col(1).isApprox(yaxis));
    BOOST_CHECK(rot.col(2).isApprox(zaxis));
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
