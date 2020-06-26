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
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include "Acts/Tests/CommonHelpers/DataDirectory.hpp"
#include "Acts/Visualization/GeometryView.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"

#include "TGeoManager.h"

namespace Acts {

namespace Test {

/// @brief struct to load the global geometry
struct RootGeometry {
  RootGeometry() {
    auto path = Acts::Test::getDataPath("panda.root");
    TGeoManager::Import(path.c_str());
  }
};

RootGeometry rGeometry = RootGeometry();

GeometryContext tgContext = GeometryContext();

/// @brief Unit test checking the match probability
BOOST_AUTO_TEST_CASE(TGeoLayerBuilderTests) {
  using TglConfig = TGeoLayerBuilder::LayerConfig;

  TglConfig b0Config;
  b0Config.layerName = "*";
  b0Config.sensorNames = {"PixelActiveo2_1", "PixelActiveo4_1",
                          "PixelActiveo5_1", "PixelActiveo6_1"};
  b0Config.localAxes = "XYZ";
  b0Config.parseRanges = {{binR, {0., 40_mm}}, {binZ, {-60_mm, 15_mm}}};
  b0Config.envelope = {0_mm, 0_mm};

  TglConfig eAllConfig;
  eAllConfig.layerName = "*";
  eAllConfig.sensorNames = {"PixelActiveo2_1", "PixelActiveo4_1",
                            "PixelActiveo5_1", "PixelActiveo6_1"};
  eAllConfig.localAxes = "XYZ";
  eAllConfig.parseRanges = {{binR, {0., 40_mm}}, {binZ, {16_mm, 60_mm}}};
  eAllConfig.splitConfigs = {{binZ, 5_mm}};
  eAllConfig.envelope = {0_mm, 0_mm};

  std::vector<TglConfig> cConfigs = {b0Config};
  std::vector<TglConfig> pConfigs = {eAllConfig};

  TGeoLayerBuilder::Config tglbConfig;
  tglbConfig.configurationName = "Pixels";
  tglbConfig.layerConfigurations[1] = cConfigs;
  tglbConfig.layerConfigurations[2] = pConfigs;

  auto surfaceArrayCreator = std::make_shared<const SurfaceArrayCreator>(
      getDefaultLogger("SurfaceArrayCreator", Logging::VERBOSE));

  LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<const LayerCreator>(
      lcConfig, getDefaultLogger("LayerCreator", Logging::VERBOSE));
  tglbConfig.layerCreator = layerCreator;

  ProtoLayerHelper::Config plhConfig;
  auto protoLayerHelper = std::make_shared<const ProtoLayerHelper>(
      plhConfig, getDefaultLogger("ProtoLayerHelper", Logging::VERBOSE));
  tglbConfig.protoLayerHelper = protoLayerHelper;

  TGeoLayerBuilder tglb(tglbConfig,
                        getDefaultLogger("TGeoLayerBuilder", Logging::VERBOSE));

  ObjVisualization objVis;

  auto centralLayers = tglb.centralLayers(tgContext);
  BOOST_TEST(centralLayers.size() = 1u);
  BOOST_TEST(tglb.detectorElements().size() == 14u);

  auto positiveLayers = tglb.positiveLayers(tgContext);
  // Check that it's split into two layers
  size_t ipl = 0;
  BOOST_TEST(positiveLayers.size() = 2u);
  BOOST_TEST(tglb.detectorElements().size() == 14u + 16u);
  for (const auto& pLayer : positiveLayers) {
    auto sArray = pLayer->surfaceArray();
    if (sArray) {
      for (auto& surface : sArray->surfaces()) {
        GeometryView::drawSurface(objVis, *surface, tgContext);
      }
    }
    objVis.write("PositiveLayer_" + std::to_string(ipl++));
    objVis.clear();
  }
}

}  // namespace Test

}  // namespace Acts