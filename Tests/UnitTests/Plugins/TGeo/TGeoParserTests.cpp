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
#include "Acts/Plugins/TGeo/TGeoParser.hpp"
#include "Acts/Plugins/TGeo/TGeoSurfaceConverter.hpp"
#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"

#include "TGeoManager.h"

namespace Acts {

namespace Test {

/// @brief struct to load the global geometry
struct RootGeometry {
  RootGeometry() {
    TGeoManager::Import("http://cern.ch/asalzbur/acts/panda.root");
  }
};

GeometryContext tgContext = GeometryContext();

RootGeometry rGeometry = RootGeometry();

/// @brief Unit test Parsing a TGeo geometry
BOOST_AUTO_TEST_CASE(TGeoParser_Pixel) {
  if (gGeoManager != nullptr) {
    std::string volumeName = "*";
    TGeoParser::Options tgpOptions;
    tgpOptions.volumeNames = {volumeName};
    tgpOptions.targetNames = {"PixelActiveo2_1", "PixelActiveo4_1",
                              "PixelActiveo5_1", "PixelActiveo6_1"};
    std::string axes = "XYZ";
    double scale = 10.;

    TGeoParser::State tgpState;
    tgpState.volume = gGeoManager->GetTopVolume();

    // Parse the full ones
    TGeoParser::select(tgpState, tgpOptions);

    // This should select 176 PixelActive modules
    BOOST_TEST(tgpState.selectedNodes.size() == 176u);

    /// Convert into surfaces using the TGeoSurfaceConverter & Draw them
    ObjVisualization objVis;
    for (auto& snode : tgpState.selectedNodes) {
      const auto& shape = *(snode.node->GetVolume()->GetShape());
      const auto& transform = *(snode.transform.get());
      auto surface =
          TGeoSurfaceConverter::toSurface(shape, transform, axes, scale);
      GeometryVisualization::drawSurface(objVis, *surface, tgContext,
                                         Transform3D::Identity(), 1, false,
                                         {120, 0, 0});
    }
    objVis.write("PixelActive");
  }
}

/// @brief Unit test Parsing a TGeo geometrys
BOOST_AUTO_TEST_CASE(TGeoParser_Pixel_SelectInnermost) {
  if (gGeoManager != nullptr) {
    std::string volumeName = "*";
    TGeoParser::Options tgpOptions;
    tgpOptions.volumeNames = {volumeName};
    tgpOptions.targetNames = {"PixelActiveo2_1", "PixelActiveo4_1",
                              "PixelActiveo5_1", "PixelActiveo6_1"};
    tgpOptions.parseRanges.push_back({binR, {0., 40.}});
    tgpOptions.parseRanges.push_back({binZ, {-60., 15.}});
    tgpOptions.unit = 10.;

    std::string axes = "XYZ";

    TGeoParser::State tgpState;
    tgpState.volume = gGeoManager->GetTopVolume();

    // Parse the full ones
    TGeoParser::select(tgpState, tgpOptions);

    // This should select 14 PixelActive modules
    BOOST_TEST(tgpState.selectedNodes.size() == 14u);

    /// Convert into surfaces using the TGeoSurfaceConverter & Draw them
    ObjVisualization objVis;
    for (auto& snode : tgpState.selectedNodes) {
      const auto& shape = *(snode.node->GetVolume()->GetShape());
      const auto& transform = *(snode.transform.get());
      auto surface = TGeoSurfaceConverter::toSurface(shape, transform, axes,
                                                     tgpOptions.unit);
      GeometryVisualization::drawSurface(objVis, *surface, tgContext,
                                         Transform3D::Identity(), 1, false,
                                         {120, 0, 0});
    }
    objVis.write("PixelActive_Innemost");
  }
}

}  // namespace Test
}  // namespace Acts
