// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"

#include <fstream>

BOOST_AUTO_TEST_SUITE(SurfaceSvgConverter)

BOOST_AUTO_TEST_CASE(PlanarSurfaces) {

  // Default Geometry context
  Acts::GeometryContext geoCtx;

  // Planar style
  Acts::Svg::Style planarStyle;
  planarStyle.fillColor = { 51, 153, 255 };
  planarStyle.highlightColor = { 255, 153, 51};
  planarStyle.highlights = { "mouseover", "mouseout" };
  planarStyle.strokeWidth = 0.5;

  // Rectangle case: Acts object
  auto rectangleBounds = std::make_shared<Acts::RectangleBounds>(36., 64.);
  auto transform = Acts::Transform3::Identity();
  transform.pretranslate(Acts::Vector3{20., 20., 100.});
  auto rectanglePlane = Acts::Surface::makeShared<Acts::PlaneSurface>(transform, rectangleBounds);
  // Svg proto object & actual object
  auto rectangleSvgT = Acts::Svg::convert(geoCtx, *rectanglePlane, 1u, planarStyle, true);
  auto rectangleXYT = Acts::Svg::xyView(rectangleSvgT, std::string("rectangleXYT"));
  Acts::Svg::toFile({rectangleXYT}, std::string("RectangleXY_asTemplate.svg"));
  // Positioned 
  auto rectangleSvg = Acts::Svg::convert(geoCtx, *rectanglePlane, 0, planarStyle, false);
  auto rectangleXY = Acts::Svg::xyView(rectangleSvg, std::string("rectangleXY"));
  Acts::Svg::toFile({rectangleXY}, std::string("RectangleXY_positioned.svg"));
  // As sheet 
  auto rectangleXYSheet = Acts::Svg::xySheet(rectangleSvgT, std::string("rectangleXYT"));
  Acts::Svg::toFile({rectangleXYSheet}, std::string("RectangleXY_sheet.svg"));






}

BOOST_AUTO_TEST_SUITE_END()
