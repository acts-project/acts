// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include <Acts/Plugins/ActSVG/LayerSvgConverter.hpp>
#include <Acts/Plugins/ActSVG/TrackingGeometrySvgConverter.hpp>

namespace {

static inline Acts::Svg::Style layerStyle() {
  Acts::Svg::Style defaultStyle;
  defaultStyle.fillColor = {51, 153, 255};
  defaultStyle.fillOpacity = 0.75;
  defaultStyle.highlightColor = {255, 153, 51};
  defaultStyle.highlights = {"mouseover", "mouseout"};
  defaultStyle.strokeColor = {25, 25, 25};
  defaultStyle.strokeWidth = 0.5;
  defaultStyle.nSegments = 72u;

  return defaultStyle;
};

static inline Acts::Svg::TrackingGeometryConverter::Options
trackingGeometryOptions() {
  Acts::GeometryIdentifier geoID(0);

  Acts::Svg::LayerConverter::Options lOptions;

  lOptions.name = "layer";
  lOptions.surfaceStyles =
      Acts::GeometryHierarchyMap<Acts::Svg::Style>({{geoID, layerStyle()}});

  Acts::Svg::TrackingGeometryConverter::Options tgOptions;
  tgOptions.prefix = "";
  tgOptions.layerOptions =
      Acts::GeometryHierarchyMap<Acts::Svg::LayerConverter::Options>(
          {{geoID, lOptions}});

  return tgOptions;
};

}  // namespace

namespace ActsExamples {

static Acts::Svg::Style s_defaultLayerStyle = layerStyle();

static Acts::Svg::TrackingGeometryConverter::Options
    s_defaultTrackingGeometryOptions = trackingGeometryOptions();
}  // namespace ActsExamples
