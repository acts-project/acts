// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/ProtoLayerHelper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/IVisualization.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"

#include <cmath>

namespace Acts {

namespace Test {
namespace Layers {
BOOST_AUTO_TEST_SUITE(Geometry)

BOOST_AUTO_TEST_CASE(ProtoLayerHelperTests) {
  ProtoLayerHelper::Config plhConfig;
  ProtoLayerHelper plHelper(plhConfig);

  GeometryContext tgContext = GeometryContext();

  ObjVisualization objVis;

  CylindricalTrackingGeometry ctGeometry(tgContext);
  CylindricalTrackingGeometry::DetectorStore dStore;

  /// Helper method to prepare the streams & helpers
  /// @param path is the file path
  /// @param clear ist he indicator to clear the helper
  auto write = [&](IVisualization& helper, const std::string& path,
                   bool clear = true) -> void {
    helper.write(path);
    if (clear) {
      helper.clear();
    }
  };

  std::vector<double> layerRadii = {32., 72., 116., 172.};
  std::vector<std::pair<int, int>> layerBinning = {
      {16, 14}, {32, 14}, {52, 14}, {78, 14}};
  std::vector<double> moduleTiltPhi = {0.145, 0.145, 0.145, 0.145};
  std::vector<double> moduleHalfX = {8.4, 8.4, 8.4, 8.4};
  std::vector<double> moduleHalfY = {36., 36., 36., 36.};
  std::vector<double> moduleThickness = {0.15, 0.15, 0.15, 0.15};

  std::vector<const Surface*> volumeSurfaces;
  for (size_t ilp = 0; ilp < layerRadii.size(); ++ilp) {
    std::vector<const Surface*> layerSurfaces = ctGeometry.surfacesCylinder(
        dStore, moduleHalfX[ilp], moduleHalfY[ilp], moduleThickness[ilp],
        moduleTiltPhi[ilp], layerRadii[ilp], 2., 5., layerBinning[ilp]);
    volumeSurfaces.insert(volumeSurfaces.begin(), layerSurfaces.begin(),
                          layerSurfaces.end());
  }

  IVisualization::ColorType unsortedColor = {252, 160, 0};
  for (auto& sf : volumeSurfaces) {
    Visualization::drawSurface(objVis, *sf, tgContext, Transform3D::Identity(),
                               1, false, unsortedColor);
  }
  // Draw the all surfaces
  write(objVis, "ProtoLayerHalper_CylinderLayers_unsorted", true);

  // Sort into ProtoLayers
  auto radialLayers = plHelper.protoLayers(tgContext, volumeSurfaces, binR, 5.);

  BOOST_CHECK(radialLayers.size() == 4);

  std::vector<IVisualization::ColorType> radiallySortedColors = {
      {102, 204, 255}, {102, 255, 153}, {255, 204, 102}, {204, 102, 0}};

  size_t il = 0;
  for (auto& layer : radialLayers) {
    for (auto& sf : layer.surfaces()) {
      Visualization::drawSurface(objVis, *sf, tgContext,
                                 Transform3D::Identity(), 1, false,
                                 radiallySortedColors[il]);
    }
    ++il;
  }

  // Draw the sorted surfaces
  write(objVis, "ProtoLayerHalper_CylinderLayers_radially", true);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Layers
}  // namespace Test

}  // namespace Acts
