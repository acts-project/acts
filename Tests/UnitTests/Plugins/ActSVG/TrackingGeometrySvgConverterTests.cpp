// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/TemporaryDirectory.hpp"
#include "ActsPlugins/ActSVG/LayerSvgConverter.hpp"
#include "ActsPlugins/ActSVG/SvgUtils.hpp"
#include "ActsPlugins/ActSVG/TrackingGeometrySvgConverter.hpp"

#include <format>
#include <fstream>
#include <memory>
#include <vector>

Acts::GeometryContext tgContext;

BOOST_AUTO_TEST_SUITE(ActSvg)

BOOST_AUTO_TEST_CASE(CylindricalTrackingGeometrySvg) {
  ActsPlugins::Svg::Style cylinderLayerStyle;
  cylinderLayerStyle.fillColor = {51, 153, 255};
  cylinderLayerStyle.fillOpacity = 0.75;
  cylinderLayerStyle.highlightColor = {255, 153, 51};
  cylinderLayerStyle.highlights = {"mouseover", "mouseout"};
  cylinderLayerStyle.strokeColor = {25, 25, 25};
  cylinderLayerStyle.strokeWidth = 0.5;
  cylinderLayerStyle.quarterSegments = 72u;

  Acts::GeometryIdentifier geoID{0};

  Acts::Test::CylindricalTrackingGeometry cGeometry(tgContext);
  auto tGeometry = cGeometry();

  ActsPlugins::Svg::LayerConverter::Options lOptions;
  lOptions.name = "cylinder_layer_";
  lOptions.surfaceStyles = Acts::GeometryHierarchyMap<ActsPlugins::Svg::Style>(
      {{geoID, cylinderLayerStyle}});

  ActsPlugins::Svg::TrackingGeometryConverter::Options tgOptions;
  tgOptions.prefix = "utest_geometry_";
  tgOptions.layerOptions =
      Acts::GeometryHierarchyMap<ActsPlugins::Svg::LayerConverter::Options>(
          {{geoID, lOptions}});

  auto geometrySheets = ActsPlugins::Svg::TrackingGeometryConverter::convert(
      tgContext, *tGeometry, tgOptions);

  for (const auto& s : geometrySheets) {
    ActsPlugins::Svg::toFile({s}, s._id + ".svg");
  }
}

BOOST_AUTO_TEST_CASE(CylindricalTrackingGeometrySvgGen3) {
  Acts::Test::CylindricalTrackingGeometry cGeometry(tgContext, true);
  Acts::Test::TemporaryDirectory tmp{};

  auto tGeometry = cGeometry();

  auto objects = ActsPlugins::Svg::drawTrackingGeometry(tgContext, *tGeometry,
                                                        actsvg::views::x_y{});

  ActsPlugins::Svg::toFile(objects, tmp.path() / "utest_geometry_gen3_xy.svg");

  objects = ActsPlugins::Svg::drawTrackingGeometry(
      tgContext, *tGeometry, actsvg::views::z_r{}, true, true);
  ActsPlugins::Svg::toFile(objects, tmp.path() / "utest_geometry_gen3_zr.svg");

  objects = ActsPlugins::Svg::drawSurfaceArrays(tgContext, *tGeometry);

  for (const auto& obj : objects) {
    ActsPlugins::Svg::toFile(
        {obj},
        tmp.path() /
            std::format("utest_geometry_gen3_surface_arrays_{}.svg", obj._id));
  }
}

BOOST_AUTO_TEST_SUITE_END()
