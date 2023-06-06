// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "ActsExamples/MultiWire/MultiWireStructureBuilder.hpp"

#include <iostream>

BOOST_AUTO_TEST_SUITE(MultiWireStructureBuilder)

BOOST_AUTO_TEST_CASE(MultiWireStructureBuilderFromGDML) {
  Acts::GeometryContext gctx;

  auto mwCfg = ActsExamples::MultiWireStructureBuilder::Config();

  mwCfg.name = "multi_layer_1";
  mwCfg.sensitiveNames = {"Inner_Skin_ml0"};
  mwCfg.passiveNames = {"xxx"};

  auto strawSurfaces = ActsExamples::MultiWireHelper::getStrawSurfaces(
      mwCfg.sensitiveNames, mwCfg.passiveNames);
  auto mwBounds =
      ActsExamples::MultiWireHelper::getMultiWireBounds(strawSurfaces);

  mwCfg.strawSurfaces = strawSurfaces;
  mwCfg.lbinning =
      ActsExamples::MultiWireHelper::layerBinning(strawSurfaces, mwBounds);
  mwCfg.multiWireBounds = mwBounds;

  Acts::Experimental::RootDetectorVolumes roots;
  auto mwBuilder = ActsExamples::MultiWireStructureBuilder(mwCfg);

  auto mlComponent = mwBuilder.construct(roots, gctx);

  BOOST_CHECK(roots.volumes.size() == 1u);
  BOOST_CHECK(mlComponent.portals.size() == 6u);

  auto edv = roots.volumes.front();
  Acts::ObjVisualization3D obj;

  Acts::GeometryView3D::drawDetectorVolume(obj, *edv, gctx);
  obj.write("1multiwire.obj");
}

BOOST_AUTO_TEST_SUITE_END()