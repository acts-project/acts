// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "ActsPlugins/ActSVG/IndexedSurfacesSvgConverter.hpp"
#include "ActsPlugins/ActSVG/SvgUtils.hpp"
#include <actsvg/meta.hpp>

#include <tuple>
#include <vector>

namespace Acts {

class SurfaceArray;

}  // namespace Acts

/// @namespace ActsPlugins::Svg::SurfaceArrayConverter
/// @ingroup actsvg_plugin
/// Utilities to convert a surface array into proto objects for SVG output

namespace ActsPlugins::Svg::SurfaceArrayConverter {

/// @ingroup actsvg_plugin
/// Nested options struct
struct Options {
  /// Hierarchy map of styles
  Acts::GeometryHierarchyMap<Style> surfaceStyles;
  /// The Grid converter options
  GridConverter::Options gridOptions;
};

/// Convert a surface array into needed constituents
/// @ingroup actsvg_plugin
///
/// @param gctx is the geometry context of the conversion call
/// @param surfaceArray is the surface to convert
/// @param cOptions the conversion options
///
/// @note the type of view is auto-generated from the binning information,
///       it transforms the surface array into an indexed array of surfaces
///       and then uses these proto objects, one can thus directly use the
///       view function of the indexed surface grid
///
/// @return a collection of proto surface object and a grid, and associations
ProtoIndexedSurfaceGrid convert(const Acts::GeometryContext& gctx,
                                const Acts::SurfaceArray& surfaceArray,
                                const Options& cOptions = Options());

}  // namespace ActsPlugins::Svg::SurfaceArrayConverter
