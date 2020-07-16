// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/TGeo/TGeoSurfaceConverter.hpp"
#include "Acts/Surfaces/PlanarBooleanBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Visualization/GeometryView.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoTrd2.h"
#include "TGeoVolume.h"
#include "TView.h"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

GeometryContext tgContext = GeometryContext();

/// @brief Unit test to convert a TGeoUnion into a Plane
///
/// * The TGeoTrd1 can only have (x/X)(z/Z) orientation
BOOST_AUTO_TEST_CASE(TGeoTrd1_to_PlaneSurface) {
  ObjVisualization objVis;

  new TGeoManager("Default", "Geometry imported from GDML");

  // TRANSFORMATION MATRICES
  // Shape: G4Trd type: TGeoTrd2
  double dx1 = 0.014249999999999999;
  double dx2 = 0.014249999999999999;
  double dy1 = 3.3076000000000003;
  double dy2 = 3.7423500000000001;
  double dz = 2.7217500000000001;
  TGeoShape *pG4TrdLeft = new TGeoTrd2("G4TrdLeft", dx1, dx2, dy1, dy2, dz);
  // Shape: G4Trd type: TGeoTrd2
  dx1 = 0.014249999999999999;
  dx2 = 0.014249999999999999;
  dy1 = 2.7867000000000002;
  dy2 = 3.3064999999999998;
  dz = 3.2542499999999999;
  TGeoShape *pG4TrdRight = new TGeoTrd2("G4TrdRight", dx1, dx2, dy1, dy2, dz);
  // Combi transformation:
  double dx = 0;
  double dy = 0;
  dz = 3.26125;
  TGeoCombiTrans *pMatrixLeft = new TGeoCombiTrans("CombiTransLeft");
  pMatrixLeft->SetTranslation(dx, dy, dz);
  // Combi transformation:
  dx = 0;
  dy = 0;
  dz = -2.7287500000000002;
  TGeoCombiTrans *pMatrixRight = new TGeoCombiTrans("CombiTransRight");
  pMatrixRight->SetTranslation(dx, dy, dz);
  TGeoUnion *pBoolNode =
      new TGeoUnion(pG4TrdLeft, pG4TrdRight, pMatrixLeft, pMatrixRight);
  // Shape: ECSensor1 type: TGeoCompositeShape
  TGeoShape *uShape = new TGeoCompositeShape("CompositeShape", pBoolNode);
  // Volume: SCT::ECSensor1
  TGeoVolume *uVolume = new TGeoVolume("Volume", uShape);
  uVolume->SetLineColor(920);
  uVolume->SetVisLeaves(kTRUE);

  gGeoManager->SetTopVolume(uVolume);

  // CLOSE GEOMETRY
  gGeoManager->CloseGeometry();

  auto uPlane = TGeoSurfaceConverter::toSurface(*uVolume->GetShape(),
                                                *gGeoIdentity, "YZX", 1);

  GeometryView::drawSurface(objVis, *uPlane, tgContext);

  objVis.write("TGeoConversion_TGeoUnion_PlaneSurface");
}

}  // namespace Test

}  // namespace Acts