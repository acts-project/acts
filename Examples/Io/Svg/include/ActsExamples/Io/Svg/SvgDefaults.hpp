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

static inline Acts::Svg::Style backgroundStyle() {
  Acts::Svg::Style bgStyle;
  bgStyle.fillColor = {55, 55, 55};
  bgStyle.fillOpacity = 0.50;
  bgStyle.highlights = {};
  bgStyle.strokeColor = {25, 25, 25};
  bgStyle.strokeWidth = 0.5;
  bgStyle.nSegments = 72u;
  return bgStyle;
};

static inline Acts::Svg::Style pointStyle() {
  Acts::Svg::Style pStyle;
  pStyle.fillColor = {200, 0, 0};
  pStyle.fillOpacity = 1.0;
  pStyle.highlightColor = {0, 200, 0};
  pStyle.highlights = {"mouseover", "mouseout"};
  pStyle.strokeColor = {0, 0, 0};
  pStyle.strokeWidth = 0.5;
  pStyle.nSegments = 72u;

  return pStyle;
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

static inline Acts::Svg::TrackingGeometryConverter::Options
backgroundGeometryOptions() {
  Acts::GeometryIdentifier geoID(0);

  Acts::Svg::LayerConverter::Options lOptions;

  lOptions.name = "layer";
  lOptions.surfaceStyles = Acts::GeometryHierarchyMap<Acts::Svg::Style>(
      {{geoID, backgroundStyle()}});
  lOptions.moduleInfo = false;
  lOptions.gridInfo = false;
  lOptions.zRange = {-100, 100};
  lOptions.phiRange = {-0.2, 0.2};

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

static Acts::Svg::Style s_pointStyle = pointStyle();

static Acts::Svg::TrackingGeometryConverter::Options
    s_defaultTrackingGeometryOptions = trackingGeometryOptions();

static Acts::Svg::TrackingGeometryConverter::Options
    s_backgroundTrackingGeometryOptions = backgroundGeometryOptions();

}  // namespace ActsExamples
