// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorInspectors/SurfaceIndexing.hpp"

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"

#include <fstream>
#include <string>

namespace ActsExamples {

std::vector<std::shared_ptr<Acts::Surface>> readSurfacesFromJson(
    const std::string& fname) {
  std::ifstream ifile(fname.c_str());
  nlohmann::json jin;
  ifile >> jin;
  // The json surface container
  auto jSurfaces = jin["entries"];
  std::vector<std::shared_ptr<Acts::Surface>> rSurfaces = {};
  rSurfaces.reserve(jSurfaces.size());
  for (const auto& js : jSurfaces) {
    rSurfaces.push_back(Acts::surfaceFromJson(js["value"]));
  }
  return rSurfaces;
}

namespace {
Acts::Svg::IndexedSurfacesConverter::Options generateDrawOptions() {
  // The converter options
  Acts::Svg::IndexedSurfacesConverter::Options isOptions;
  // Sensitive surface stuyle
  Acts::Svg::Style sensitiveStyle;
  sensitiveStyle.fillColor = {51, 153, 255};
  sensitiveStyle.fillOpacity = 0.9;
  sensitiveStyle.highlightColor = {255, 153, 51};
  sensitiveStyle.highlights = {"onmouseover", "onmouseout"};
  sensitiveStyle.strokeWidth = 0.5;
  sensitiveStyle.strokeColor = {0, 0, 0};
  sensitiveStyle.nSegments = 72u;
  std::pair<Acts::GeometryIdentifier, Acts::Svg::Style> allSensitives = {
      Acts::GeometryIdentifier(0u), sensitiveStyle};

  // Hierarchy map of styles
  Acts::GeometryHierarchyMap<Acts::Svg::Style> surfaceStyles({allSensitives});
  isOptions.surfaceStyles = surfaceStyles;

  // The grid style
  Acts::Svg::GridConverter::Options gridOptions;
  Acts::Svg::Style gridStyle;
  gridStyle.fillOpacity = 0.;
  gridStyle.strokeColor = {0, 0, 255};
  gridStyle.strokeWidth = 1.;
  gridStyle.highlightStrokeWidth = 3;
  gridStyle.highlightStrokeColor = {255, 0, 0};
  gridOptions.style = gridStyle;

  isOptions.gridOptions = gridOptions;
  return isOptions;
};

}  // namespace

// Cylindrical inspector
CylindricalDetectorIndexing::CylindricalDetectorIndexing(
    const std::string& fname)
    : SurfaceIndexing<2u>(fname) {}

// Run the inspector code
void CylindricalDetectorIndexing::inspect(
    const std::string& name,
    const std::array<std::array<Acts::ActsScalar, 2u>, 2u> qRange,
    const std::vector<SurfaceIndexing<2u>::Binning>& binnings,
    const std::vector<SurfaceIndexing<2u>::Support>& supports) {
  SurfaceIndexing<2u>::createIndexing(name, qRange, binnings, supports,
                                      generateDrawOptions());
}

}  // namespace ActsExamples
