// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include "actsvg/core.hpp"
#include "actsvg/display/geometry.hpp"

namespace Acts {

using Point3Container = std::vector<Vector3>;

class Surface;

namespace Svg {

using ProtoSurface = actsvg::proto::surface<Point3Container>;

/// Convert into a svg::proto surface
///
/// @param gtcx is the geometry context of the conversion call
/// @param surface is the surface to convert
/// @param surfaceStyle is the surface style in question
/// @param templateSurface is a directive to create a template surface w/o transform
///
/// @return a proto surface object
ProtoSurface convert(const GeometryContext& gctx, const Surface& surface,
                     const Style& surfaceStyle,
                     bool templateSurface = false);

/// Convert into an acts::svg::object with an XY view
///
/// @param pSurface is the proto object
/// @param identification is the to be translated id_ for actsvg
///
/// @return an svg object that can be written out directly to disc
static inline actsvg::svg::object surfaceViewXY(const ProtoSurface& pSurface,
                                         const std::string& identification) {
  actsvg::views::x_y xyView;
  return actsvg::display::surface(identification, pSurface, xyView, true);
}

/// Convert into an acts::svg::object with an XY sheet
///
/// @param pSurface is the proto object
/// @param identification is the to be translated id_ for actsvg
///
/// @return an svg object that can be written out directly to disc
static inline actsvg::svg::object surfaceSheetXY(const ProtoSurface& pSurface,
                                         const std::string& identification) {
  return actsvg::display::surface_sheet_xy(identification, pSurface);
}

}  // namespace Svg

}  // namespace Acts
