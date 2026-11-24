// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Root/BlueprintBuilder.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <stdexcept>
#include <string>

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

static_assert(detail::HasPassiveVolume<TGeoBlueprintBuilderBackend>);

namespace {

/// Builds a small TGeo world (in cm) with a support tube, an endplate disc
/// and a box.
void buildPassiveGeometry() {
  new TGeoManager("passive", "passive elements");
  auto* mat = new TGeoMaterial("Al", 26.98, 13, 2.7);
  auto* med = new TGeoMedium("MED", 1, mat);
  TGeoVolume* top = gGeoManager->MakeBox("World", med, 500., 500., 500.);
  gGeoManager->SetTopVolume(top);

  // Tube-like: radial thickness 0.4 cm, z half-length 100 cm.
  TGeoVolume* tube = gGeoManager->MakeTube("SupportTube", med, 20., 20.4, 100.);
  top->AddNode(tube, 1);

  // Disc-like: radial extent 10 cm, z half-length 0.5 cm, shifted to z=150 cm.
  TGeoVolume* disc = gGeoManager->MakeTube("Endplate", med, 5., 15., 0.5);
  top->AddNode(disc, 1, new TGeoTranslation(0., 0., 150.));

  // Neither tube- nor disc-like.
  TGeoVolume* box = gGeoManager->MakeBox("Box", med, 1., 1., 1.);
  top->AddNode(box, 1);

  gGeoManager->CloseGeometry();
}

}  // namespace

BOOST_AUTO_TEST_SUITE(RootSuite)

BOOST_AUTO_TEST_CASE(TGeoBlueprintBuilderPassiveElement) {
  buildPassiveGeometry();

  TGeoBlueprintBuilderBackend::Config cfg;
  cfg.root = gGeoManager->GetTopNode();
  cfg.lengthScale = 10.;  // cm -> mm
  ActsPlugins::BlueprintBuilder builder{cfg};

  auto gctx = GeometryContext::dangerouslyDefaultConstruct();

  // Tube-like element becomes a cylinder shell around the beam axis.
  auto tubeElement = builder.findDetElementByName("SupportTube");
  BOOST_REQUIRE(tubeElement.has_value());
  auto tubeNode = builder.makePassiveElement(*tubeElement);
  BOOST_CHECK_EQUAL(tubeNode->name(), "SupportTube");
  const Volume& tubeVolume = tubeNode->build(BlueprintOptions{}, gctx);
  const auto* tubeBounds =
      dynamic_cast<const CylinderVolumeBounds*>(&tubeVolume.volumeBounds());
  BOOST_REQUIRE(tubeBounds != nullptr);
  CHECK_CLOSE_ABS(tubeBounds->get(CylinderVolumeBounds::eMinR), 200., 1e-9);
  CHECK_CLOSE_ABS(tubeBounds->get(CylinderVolumeBounds::eMaxR), 204., 1e-9);
  CHECK_CLOSE_ABS(tubeBounds->get(CylinderVolumeBounds::eHalfLengthZ), 1000.,
                  1e-9);
  CHECK_CLOSE_ABS(tubeVolume.localToGlobalTransform(gctx).translation().z(), 0.,
                  1e-9);

  // Disc-like element becomes an annular disc at its placement, with an
  // explicit volume name.
  auto discElement = builder.findDetElementByName("Endplate");
  BOOST_REQUIRE(discElement.has_value());
  auto discNode = builder.makePassiveElement(*discElement, "PositiveEndplate");
  BOOST_CHECK_EQUAL(discNode->name(), "PositiveEndplate");
  const Volume& discVolume = discNode->build(BlueprintOptions{}, gctx);
  const auto* discBounds =
      dynamic_cast<const CylinderVolumeBounds*>(&discVolume.volumeBounds());
  BOOST_REQUIRE(discBounds != nullptr);
  CHECK_CLOSE_ABS(discBounds->get(CylinderVolumeBounds::eMinR), 50., 1e-9);
  CHECK_CLOSE_ABS(discBounds->get(CylinderVolumeBounds::eMaxR), 150., 1e-9);
  CHECK_CLOSE_ABS(discBounds->get(CylinderVolumeBounds::eHalfLengthZ), 5.,
                  1e-9);
  CHECK_CLOSE_ABS(discVolume.localToGlobalTransform(gctx).translation().z(),
                  1500., 1e-9);

  // A shape that is neither tube- nor disc-like is rejected.
  auto boxElement = builder.findDetElementByName("Box");
  BOOST_REQUIRE(boxElement.has_value());
  BOOST_CHECK_THROW(builder.makePassiveElement(*boxElement),
                    std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
