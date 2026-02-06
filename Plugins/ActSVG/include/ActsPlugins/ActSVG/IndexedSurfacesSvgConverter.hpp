// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/IndexGrid.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "Acts/Utilities/TypeList.hpp"
#include "ActsPlugins/ActSVG/GridSvgConverter.hpp"
#include "ActsPlugins/ActSVG/SurfaceSvgConverter.hpp"
#include "ActsPlugins/ActSVG/SvgUtils.hpp"
#include <actsvg/core.hpp>
#include <actsvg/meta.hpp>

#include <tuple>
#include <vector>

namespace ActsPlugins::Svg {

using ProtoSurface = actsvg::proto::surface<std::vector<Acts::Vector3>>;
using ProtoSurfaces = std::vector<ProtoSurface>;
using ProtoGrid = actsvg::proto::grid;
using ProtoAssociations = std::vector<std::vector<std::size_t>>;

using ProtoIndexedSurfaceGrid =
    std::tuple<ProtoSurfaces, ProtoGrid, ProtoAssociations>;

/// @ingroup actsvg_plugin
namespace IndexedSurfacesConverter {
/// @ingroup actsvg_plugin
/// Nested options struct
struct Options {
  /// Hierarchy map of styles
  Acts::GeometryHierarchyMap<Style> surfaceStyles;
  /// The Grid converter options
  GridConverter::Options gridOptions;
};

/// Convert a surface array into needed constituents
///
/// @note actual conversion implementation, bottom of unrolling loop
///
/// @param gctx is the geometry context of the conversion call
/// @param surfaces the container of surfaces
/// @param indexGrid the indexGrid delegate
/// @param cOptions the conversion options
///
/// @return a collection of proto surface object and a grid, and associations
template <typename surface_container, typename index_grid>
ProtoIndexedSurfaceGrid convertImpl(const Acts::GeometryContext& gctx,
                                    const surface_container& surfaces,
                                    const index_grid& indexGrid,
                                    const Options& cOptions) {
  // The return surfaces
  ProtoSurfaces pSurfaces;
  pSurfaces.reserve(surfaces.size());

  // Make a local copy of the grid Options, in case we have to estimate
  // an additional bound
  GridConverter::Options gridOptions = cOptions.gridOptions;
  Acts::Extent constrain;

  // Estimate the radial extension
  // - for 1D phi
  // - for 2D z-phi or phi-z
  bool estimateR = (index_grid::grid_type::DIM == 1 &&
                    indexGrid.casts[0u] == Acts::AxisDirection::AxisPhi) ||
                   (index_grid::grid_type::DIM == 2 &&
                    (indexGrid.casts[0u] == Acts::AxisDirection::AxisPhi ||
                     indexGrid.casts[1u] == Acts::AxisDirection::AxisPhi));

  for (auto [is, s] : Acts::enumerate(surfaces)) {
    // Create the surface converter options
    SurfaceConverter::Options sOptions;
    // Try to see if you have predefined surface styles
    auto sfIter = cOptions.surfaceStyles.find(s->geometryId());
    if (sfIter != cOptions.surfaceStyles.end()) {
      sOptions.style = *sfIter;
    }
    auto pSurface = SurfaceConverter::convert(gctx, *s, sOptions);
    // Run the r estimation
    if (estimateR) {
      auto sExtent = s->polyhedronRepresentation(gctx, 4u).extent();
      if constexpr (index_grid::grid_type::DIM == 2u) {
        pSurface._radii[0u] =
            static_cast<float>(sExtent.medium(Acts::AxisDirection::AxisR));
      }
      constrain.extend(sExtent, {Acts::AxisDirection::AxisR});
    }
    // Add center info
    std::string centerInfo = " - center = (";
    const Acts::Vector3& center = s->center(gctx);
    centerInfo +=
        std::to_string(Acts::VectorHelpers::cast(center, indexGrid.casts[0u]));
    if (indexGrid.casts.size() > 1u) {
      centerInfo += ", ";
      centerInfo += std::to_string(
          Acts::VectorHelpers::cast(center, indexGrid.casts[1u]));
      centerInfo += ")";
    }
    pSurface._aux_info["center"] = {centerInfo};
    // Add the center info
    pSurfaces.push_back(pSurface);
  }

  // Adjust the grid options
  if constexpr (index_grid::grid_type::DIM == 1u) {
    if (indexGrid.casts[0u] == Acts::AxisDirection::AxisPhi) {
      auto estRangeR = constrain.range(Acts::AxisDirection::AxisR);
      std::array<double, 2u> rRange = {estRangeR.min(), estRangeR.max()};
      gridOptions.optionalBound = {rRange, Acts::AxisDirection::AxisR};
    }
  }

  // Create the grid
  ProtoGrid pGrid =
      GridConverter::convert(indexGrid.grid, indexGrid.casts, gridOptions);

  auto axes = indexGrid.grid.axes();

  // Specify the highlight indices
  ProtoAssociations highlightIndices;

  // 1D connections
  if constexpr (index_grid::grid_type::DIM == 1u) {
    const auto& binEdges = axes[0u]->getBinEdges();
    for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      typename index_grid::grid_type::index_t lbin;
      lbin[0u] = ib0;
      highlightIndices.push_back(indexGrid.grid.atLocalBins(lbin));
      // Register the bin naming
      std::string binInfo =
          std::string("- bin : [") + std::to_string(ib0) + std::string("]");
      double binCenter = 0.5 * (binEdges[ib0] + binEdges[ib0 - 1u]);
      binInfo += "\n - center : (" + std::to_string(binCenter) + ")";
      pGrid._bin_ids.push_back(binInfo);
    }
  }
  // 2D connections
  if constexpr (index_grid::grid_type::DIM == 2u) {
    const auto& binEdges0 = axes[0u]->getBinEdges();
    const auto& binEdges1 = axes[1u]->getBinEdges();
    for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      for (unsigned int ib1 = 1u; ib1 <= axes[1u]->getNBins(); ++ib1) {
        typename index_grid::grid_type::index_t lbin;
        lbin[0u] = ib0;
        lbin[1u] = ib1;
        highlightIndices.push_back(indexGrid.grid.atLocalBins(lbin));
        // Register the bin naming
        std::string binInfo = std::string("- bin : [") + std::to_string(ib0) +
                              std::string(", ") + std::to_string(ib1) +
                              std::string("]");
        double binCenter0 = 0.5 * (binEdges0[ib0] + binEdges0[ib0 - 1u]);
        double binCenter1 = 0.5 * (binEdges1[ib1] + binEdges1[ib1 - 1u]);
        binInfo += "\n - center : (" + std::to_string(binCenter0) + ", " +
                   std::to_string(binCenter1) + ")";
        pGrid._bin_ids.push_back(binInfo);
        if (estimateR) {
          pGrid._reference_r =
              static_cast<float>(constrain.medium(Acts::AxisDirection::AxisR));
        }
      }
    }
  }
  return {pSurfaces, pGrid, highlightIndices};
}

}  // namespace IndexedSurfacesConverter

/// @ingroup actsvg_plugin
namespace View {

/// Convert into an ActsPlugins::Svg::object with an XY view
///
/// @param pIndexGrid is the proto object
/// @param identification is the to be translated id_ for actsvg
///
/// @return an svg object that can be written out directly
static inline actsvg::svg::object xy(const ProtoIndexedSurfaceGrid& pIndexGrid,
                                     const std::string& identification) {
  auto [pSurfaces, pGrid, pIndices] = pIndexGrid;

  // Create the master object
  actsvg::svg::object xyIndexedGrid;
  xyIndexedGrid._id = identification;
  xyIndexedGrid._tag = "g";

  // The surfaces
  std::vector<actsvg::svg::object> sObs;
  for (const auto& s : pSurfaces) {
    if (pGrid._type == ProtoGrid::e_z_phi) {
      sObs.push_back(zrphi(s, s._name));
    } else {
      sObs.push_back(xy(s, s._name));
    }
  }

  // The grid as provided
  auto gOb =
      actsvg::display::grid(identification + std::string("_grid"), pGrid);

  // The connectors
  actsvg::connectors::connect_action(gOb._sub_objects, sObs, pIndices);

  // Add them all
  xyIndexedGrid.add_objects(sObs);

  auto xmax = xyIndexedGrid._x_range[1u];
  // The association info boxes
  for (auto [ig, gTile] : Acts::enumerate(gOb._sub_objects)) {
    // Target surface text
    std::vector<std::string> binText;
    binText.push_back("Source:");
    binText.push_back(pGrid._bin_ids[ig]);
    binText.push_back("Target:");
    for (const auto [is, sis] : Acts::enumerate(pIndices[ig])) {
      const auto& ps = pSurfaces[sis];
      std::string oInfo = std::string("- object: ") + std::to_string(sis);
      if (ps._aux_info.contains("center")) {
        for (const auto& ci : ps._aux_info.at("center")) {
          oInfo += ci;
        }
      }
      binText.push_back(oInfo);
    }
    // Make the connected text
    std::string cTextId =
        identification + std::string("_ct_") + std::to_string(ig);
    auto cText = actsvg::draw::connected_text(
        cTextId, {static_cast<actsvg::scalar>(1.1 * xmax), 0}, binText,
        actsvg::style::font{}, actsvg::style::transform{}, gTile);
    xyIndexedGrid.add_object(cText);
  }
  xyIndexedGrid.add_object(gOb);
  // Return the grid
  return xyIndexedGrid;
}

/// Fake zphi method, delegate to xy
///
/// @param pIndexGrid is the proto object
/// @param identification is the to be translated id_ for actsvg
///
/// @note this works, because the actual display at the end is 2D
/// in an x-y plane and the parameters are appropriately defined
///
/// @return an svg object that can be written out directly
static inline actsvg::svg::object zphi(
    const ProtoIndexedSurfaceGrid& pIndexGrid,
    const std::string& identification) {
  return xy(pIndexGrid, identification);
}

}  // namespace View
}  // namespace ActsPlugins::Svg
