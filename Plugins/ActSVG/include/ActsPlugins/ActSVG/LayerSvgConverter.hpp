// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "ActsPlugins/ActSVG/SvgUtils.hpp"
#include <actsvg/core.hpp>
#include <actsvg/meta.hpp>

#include <numbers>

namespace Acts {

class Layer;

};

namespace ActsPlugins::Svg {

using ProtoVolume = actsvg::proto::volume<std::vector<Acts::Vector3>>;

static const std::array<double, 2> noLimitZ = {
    std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()};

static const std::array<double, 2> noLimitPhi = {-std::numbers::pi,
                                                 std::numbers::pi};

/// @ingroup actsvg_plugin
namespace LayerConverter {

/// @ingroup actsvg_plugin
/// The enumeration for sheets
enum Sheets {
  eModuleInfo = 0,
  eGridInfo = 1,
  eCrossSectionXY = 2,
  eCrossSectionZR = 3
};

/// @ingroup actsvg_plugin
/// A nested options class for the layer conversion
struct Options {
  /// The name for the conversion object
  std::string name = "";
  /// The style of the surface objects
  Acts::GeometryHierarchyMap<Style> surfaceStyles;
  /// The z limit for projections
  std::array<double, 2> zRange = noLimitZ;
  /// The phi limit for projections
  std::array<double, 2> phiRange = noLimitPhi;
  /// Configuration of the views
  /// @deprecated actsvg 0.4.57 removed the grid sheet display. Kept so the
  /// option still compiles; convert() throws if it is switched on. No
  /// [[deprecated]] attribute: on a data member it also fires on every copy
  /// of Options, which says nothing about whether the option is used.
  bool gridInfo = false;
  /// Include module information
  /// @deprecated actsvg 0.4.57 removed the module sheet display; see gridInfo.
  bool moduleInfo = false;
  /// Include projection information
  bool projectionInfo = true;
  /// Label checks
  bool labelProjection = false;
  /// Label gauge for display
  double labelGauge = 0.;
};

/// Write/create the layer sheets for a given layer
///
/// @param gctx the geometry context
/// @param layer the layer to be displayed
/// @param cOptions the conversion objects
///
/// @return a vector of svg objects
std::vector<actsvg::svg::object> convert(const Acts::GeometryContext& gctx,
                                         const Acts::Layer& layer,
                                         const Options& cOptions);

}  // namespace LayerConverter

}  // namespace ActsPlugins::Svg
