// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"

namespace Acts {

class Surface;

namespace Svg {

using ProtoGrid = actsvg::proto::grid;

namespace GridConverter {

/// Nested Options struct
struct Options {
  /// A The style for the surfaces
  Style style{{150, 150, 150}, 0.5};
  /// ACTS log level
  Logging::Level logLevel = Logging::INFO;
};

/// Convert an ACTS grid into a actsvg protogrid, it currently works with
///
/// - 1D: [ binX ] , [ binY ], [ binR ] , [ binPhi ]
/// - 2D: [ binX, binY ], [ binZ, binPhi ], [ binR, binPhi ]
///
/// @tparam grid_type is the type of the grid to be converted
///
/// @param grid the grid to be converted
/// @param bValues the binning values identifying the axes
/// @param cOptions the conversion options
///
/// @return an ACTSVG proto grid for displaying
template <typename grid_type>
ProtoGrid convert(const grid_type& grid,
                  const std::array<BinningValue, grid_type::DIM>& bValues,
                  const GridConverter::Options& cOptions) {
  // The return object
  ProtoGrid pGrid;

  // Grid axes
  auto axes = grid.axes();

  /// Helper method to convert from ACTS to Grid edges
  ///
  /// @param acrtsEdges the grid edges from ACTS to be converted via static_cast
  auto convertGridEdges = [](const std::vector<Acts::ActsScalar>& actsEdges)
      -> std::vector<actsvg::scalar> {
    std::vector<actsvg::scalar> svgEdges;
    svgEdges.reserve(actsEdges.size());
    for (const auto ae : actsEdges) {
      svgEdges.push_back(static_cast<actsvg::scalar>(ae));
    }
    return svgEdges;
  };

  // The endges values
  std::vector<Acts::ActsScalar> edges0;
  std::vector<Acts::ActsScalar> edges1;

  // 1D
  if constexpr (grid_type::DIM == 1u) {
    if (bValues[0u] == binPhi and
        axes[0]->getBoundaryType() == detail::AxisBoundaryType::Closed) {
      edges1 = axes[0]->getBinEdges();
      pGrid._type = actsvg::proto::grid::e_r_phi;
    }
  }

  // 2D cases
  if constexpr (grid_type::DIM == 2u) {
    // Walk throuth the binning and translate
    if (bValues[0] == binPhi and bValues[1] == binZ) {
      //  flip to fit with actsvg convention
      edges1 = axes[0]->getBinEdges();
      edges0 = axes[1]->getBinEdges();
      pGrid._type = actsvg::proto::grid::e_z_phi;
    } else if (bValues[0] == binPhi and bValues[1] == binR) {
      //  flip to fit with actsvg convention
      edges1 = axes[0]->getBinEdges();
      edges0 = axes[1]->getBinEdges();
      pGrid._type = actsvg::proto::grid::e_r_phi;
    } else if (bValues[0] == binZ and bValues[1] == binPhi) {
      // good
      edges0 = axes[0]->getBinEdges();
      edges1 = axes[1]->getBinEdges();
      pGrid._type = actsvg::proto::grid::e_z_phi;
    } else if (bValues[0] == binR and bValues[1] == binPhi) {
      // good
      edges0 = axes[0]->getBinEdges();
      edges1 = axes[1]->getBinEdges();
      pGrid._type = actsvg::proto::grid::e_r_phi;
    } else if (bValues[0] == binX and bValues[1] == binY) {
      // good
      edges0 = axes[0]->getBinEdges();
      edges1 = axes[1]->getBinEdges();
      pGrid._type = actsvg::proto::grid::e_x_y;
    }
  }

  // Assign grid edges
  pGrid._edges_0 = convertGridEdges(edges0);
  pGrid._edges_1 = convertGridEdges(edges1);

  auto [fill, stroke] = cOptions.style.fillAndStroke();
  pGrid._fill = fill;
  pGrid._stroke = stroke;

  return pGrid;
}

}  // namespace GridConverter
}  // namespace Svg
}  // namespace Acts
