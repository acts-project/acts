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
#include "actsvg/meta.hpp"

#include <tuple>
#include <vector>

namespace Acts {

class SurfaceArray;

namespace Svg {

using ProtoSurface = actsvg::proto::surface<std::vector<Vector3>>;
using ProtoSurfaces = std::vector<ProtoSurface>;
using ProtoGrid = actsvg::proto::grid;

/// Convert a surface array into needed constituents
///
/// @param gtcx is the geometry context of the conversion call
/// @param surfaceArray is the surface to convert
/// @param surfaceStyle is the surface style in question
///
/// @note the type of view is auto-generated from the binning information
///
/// @return a collection of proto surface object and a grid with associations
std::tuple<ProtoSurfaces, ProtoGrid> convert(const GeometryContext& gctx,
                                             const SurfaceArray& surfaceArray,
                                             const Style& surfaceStyle);

}  // namespace Svg

}  // namespace Acts
