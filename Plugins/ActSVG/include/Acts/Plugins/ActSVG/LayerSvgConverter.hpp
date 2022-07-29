// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"

namespace Acts {

class Layer;

namespace Svg {

using ProtoVolume = actsvg::proto::volume<std::vector<Vector3>>;

/// Write/create the layer sheets for a given layer
///
/// @param gctx the geometry context
/// @param layer the layer to be displayed
/// @param surfaceStyle the style of the surface drawings
/// @param name an optional name of the objects
///
/// @return a vector of svg objects
std::vector<actsvg::svg::object> layerSheets(const GeometryContext& gctx,
                                             const Layer& layer,
                                             const std::string& layerName,
                                             const Style& surfaceStyle);

}  // namespace Svg

}  // namespace Acts
