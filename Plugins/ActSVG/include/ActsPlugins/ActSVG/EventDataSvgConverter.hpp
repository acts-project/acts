// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/ActSVG/SvgUtils.hpp"
#include <actsvg/core.hpp>

/// @namespace ActsPlugins::Svg::EventDataConverter
/// @ingroup actsvg_plugin
/// Utilities to convert event data to SVG objects

namespace ActsPlugins::Svg::EventDataConverter {

/// @addtogroup actsvg_plugin
/// @{

/// Write/create a 3D point in XY view
///
/// @param pos the position
/// @param size the size of the object
/// @param style the style of the object
/// @param idx the running index
///
/// @return a vector of svg objects
actsvg::svg::object pointXY(const Acts::Vector3& pos, double size,
                            const Style& style, unsigned int idx = 0);

/// Write/create a 3D point in ZR view
///
/// @param pos the position
/// @param size the size of the object
/// @param style the style of the object
/// @param idx the running index
///
/// @return a vector of svg objects
actsvg::svg::object pointZR(const Acts::Vector3& pos, double size,
                            const Style& style, unsigned int idx = 0);

/// Write/create a 3D point in a given view
///
/// @param pos the position
/// @param size the size of the object
/// @param style the style of the object
/// @param idx the running index
///
/// @return a vector of svg objects
template <typename view_type>
actsvg::svg::object point(const Acts::Vector3& pos, double size,
                          const Style& style, unsigned int idx) {
  view_type view;
  std::vector<Acts::Vector3> ps = {pos};
  auto ppos = view(ps)[0];
  auto [fill, stroke] = style.fillAndStroke();
  auto circle =
      actsvg::draw::circle("p_" + std::to_string(idx), ppos,
                           static_cast<actsvg::scalar>(size), fill, stroke);
  return circle;
}

/// @}

}  // namespace ActsPlugins::Svg::EventDataConverter
