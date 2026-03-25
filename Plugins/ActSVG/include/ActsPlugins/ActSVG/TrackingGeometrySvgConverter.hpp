// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsPlugins/ActSVG/LayerSvgConverter.hpp"
#include "ActsPlugins/ActSVG/SvgUtils.hpp"
#include <actsvg/core.hpp>
#include <actsvg/meta.hpp>

namespace Acts {

class TrackingGeometry;
class TrackingVolume;

}  // namespace Acts

namespace ActsPlugins::Svg {

using ProtoPortal = actsvg::proto::portal<std::vector<Acts::Vector3>>;
using ProtoLink = ProtoPortal::link;

/// @ingroup actsvg_plugin
namespace TrackingGeometryConverter {

/// @ingroup actsvg_plugin
/// Nested Options struct for the writing configuration
struct Options {
  /// Prefix the output names
  std::string prefix = "";
  /// Write the layer conversion options
  Acts::GeometryHierarchyMap<LayerConverter::Options> layerOptions;
};

/// @ingroup actsvg_plugin
/// State object to collect geometry-wise information
struct State {
  /// XY cross section views
  std::vector<actsvg::svg::object> xyCrossSection;
  /// ZR cross section views
  std::vector<actsvg::svg::object> zrCrossSection;

  /// Final rendered views
  std::vector<actsvg::svg::object> finalViews;
};

/// Write svg files per tracking geometry
///
/// @param gctx the geometry context
/// @param tGeometry the tracking geometry
/// @param cOptions the conversion options
///
/// @return a vector of svg objects
std::vector<actsvg::svg::object> convert(
    const Acts::GeometryContext& gctx, const Acts::TrackingGeometry& tGeometry,
    const Options& cOptions);

/// Recursivele conversion of volumes
///
/// @param gctx the geometry context
/// @param tVolume the tracking volume
/// @param cOptions the conversion options
/// @param cState [in,out] the conversion state collecting the input
///
void convert(const Acts::GeometryContext& gctx,
             const Acts::TrackingVolume& tVolume, const Options& cOptions,
             State& cState);

}  // namespace TrackingGeometryConverter

/// @ingroup actsvg_plugin
namespace TrackingGeometryProjections {

/// Options for tracking geometry projections
struct Options {
  /// Prefix for output names
  std::string prefix = "";

  /// Tracking geometry converter options
  TrackingGeometryConverter::Options trackingGeometryOptions;

  /// RZ axes ranges
  std::array<std::array<double, 2>, 2> rzAxes{};
  /// RZ eta line positions
  std::vector<double> rzEtaLines;
};

/// Convert into xy and zr projections only
///
/// @param gctx the geometry context
/// @param tGeometry the tracking volume
/// @param cOptions the conversion options
///
/// @note best performant if configured with options
/// that omit the module info and grid info
///
std::array<actsvg::svg::object, 2> convert(
    const Acts::GeometryContext& gctx, const Acts::TrackingGeometry& tGeometry,
    const Options& cOptions);
}  // namespace TrackingGeometryProjections

[[nodiscard("Not drawing svg outputs")]]
std::vector<actsvg::svg::object> drawTrackingGeometry(
    const Acts::GeometryContext& gctx, const Acts::TrackingGeometry& tGeometry,
    std::variant<actsvg::views::x_y, actsvg::views::z_r> view,
    bool drawSurfaces = true, bool highlightMaterial = false);

[[nodiscard("Not drawing svg outputs")]]
std::vector<actsvg::svg::object> drawSurfaceArrays(
    const Acts::GeometryContext& gctx, const Acts::TrackingGeometry& tGeometry);

}  // namespace ActsPlugins::Svg
