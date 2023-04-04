// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/GridAxisGenerators.hpp"
#include "Acts/Detector/IndexedSurfacesGenerator.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Plugins/ActSVG/GridSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/detail/Grid.hpp"
#include "actsvg/meta.hpp"

#include <tuple>
#include <vector>

namespace Acts {

namespace Experimental {

using namespace GridAxisGenerators;

// Helper tuple for unfolding
//
// This is a bit boiler plate, and one could make the compiler create those
// but this would increase compile time significantly
inline static auto createAxisTuple() {
  // All 1D Eq options
  EqBound eb{{0, 1}, 1};
  EqOpen eo{{0, 1}, 1};
  EqClosed ec{{0, 1}, 1};

  // All 1D Var  options
  VarBound vb{{0., 1.}};
  VarOpen vo{{0., 1.}};
  VarClosed vc{{0., 1.}};

  // All 2D EqEq options
  EqBoundEqBound ebeb{{0., 1.}, 1u, {0., 1.}, 1u};
  EqBoundEqOpen ebeo{{0., 1.}, 1u, {0., 1.}, 1u};
  EqBoundEqClosed ebec{{0., 1.}, 1u, {0., 1.}, 1u};
  EqOpenEqBound eoeb{{0., 1.}, 1u, {0., 1.}, 1u};
  EqOpenEqOpen eoeo{{0., 1.}, 1u, {0., 1.}, 1u};
  EqOpenEqClosed eoec{{0., 1.}, 1u, {0., 1.}, 1u};
  EqClosedEqBound eceb{{0., 1.}, 1u, {0., 1.}, 1u};
  EqClosedEqOpen eceo{{0., 1.}, 1u, {0., 1.}, 1u};
  EqClosedEqClosed ecec{{0., 1.}, 1u, {0., 1.}, 1u};

  // All 2D EqVar options
  EqBoundVarBound ebvb{{0., 1.}, 1u, {0., 1.}};
  EqBoundVarOpen ebvo{{0., 1.}, 1u, {0., 1.}};
  EqBoundVarClosed ebvc{{0., 1.}, 1u, {0., 1.}};
  EqOpenVarBound eovb{{0., 1.}, 1u, {0., 1.}};
  EqOpenVarOpen eovo{{0., 1.}, 1u, {0., 1.}};
  EqOpenVarClosed eovc{{0., 1.}, 1u, {0., 1.}};
  EqClosedVarBound ecvb{{0., 1.}, 1u, {0., 1.}};
  EqClosedVarOpen ecvo{{0., 1.}, 1u, {0., 1.}};
  EqClosedVarClosed ecvc{{0., 1.}, 1u, {0., 1.}};

  // All 2D VarEq options
  VarBoundEqBound vbeb{{0., 1.}, {0., 1.}, 1u};
  VarBoundEqOpen vbeo{{0., 1.}, {0., 1.}, 1u};
  VarBoundEqClosed vbec{{0., 1.}, {0., 1.}, 1u};
  VarOpenEqBound voeb{{0., 1.}, {0., 1.}, 1u};
  VarOpenEqOpen voeo{{0., 1.}, {0., 1.}, 1u};
  VarOpenEqClosed voec{{0., 1.}, {0., 1.}, 1u};
  VarClosedEqBound vceb{{0., 1.}, {0., 1.}, 1u};
  VarClosedEqOpen vceo{{0., 1.}, {0., 1.}, 1u};
  VarClosedEqClosed vcec{{0., 1.}, {0., 1.}, 1u};

  // All 2D VarEq options
  VarBoundVarBound vbvb{{0., 1.}, {0., 1.}};
  VarBoundVarOpen vbvo{{0., 1.}, {0., 1.}};
  VarBoundVarClosed vbvc{{0., 1.}, {0., 1.}};
  VarOpenVarBound vovb{{0., 1.}, {0., 1.}};
  VarOpenVarOpen vovo{{0., 1.}, {0., 1.}};
  VarOpenVarClosed vovc{{0., 1.}, {0., 1.}};
  VarClosedVarBound vcvb{{0., 1.}, {0., 1.}};
  VarClosedVarOpen vcvo{{0., 1.}, {0., 1.}};
  VarClosedVarClosed vcvc{{0., 1.}, {0., 1.}};

  return std::tie(eb, eo, ec, vb, vo, vc, ebeb, ebeo, ebec, eoeb, eoeo, eoec,
                  eceb, eceo, ecec, ebvb, ebvo, ebvc, eovb, eovo, eovc, ecvb,
                  ecvo, ecvc, vbeb, vbeo, vbec, voeb, voeo, voec, vceb, vceo,
                  vcec, vbvb, vbvo, vbvc, vovb, vovo, vovc, vcvb, vcvo, vcvc);
}

static auto s_possibleAxes = createAxisTuple();

}  // namespace Experimental

namespace Svg {

using ProtoSurface = actsvg::proto::surface<std::vector<Vector3>>;
using ProtoGrid = actsvg::proto::grid;
using ProtoIndexedSurfaceGrid =
    std::tuple<std::vector<ProtoSurface>, ProtoGrid,
               std::vector<std::vector<std::size_t>>>;

namespace IndexedSurfacesConverter {
/// Nested options struct
struct Options {
  /// Hierarchy map of styles
  GeometryHierarchyMap<Style> surfaceStyles;

  /// The Grid converter options
  GridConverter::Options gridOptions;

  /// ACTS Logging level
  Logging::Level logLevel = Logging::INFO;
};

/// @brief Convert the single delegate if it is of the type of the reference
///
/// @note It will do nothing if the type does not match
///
/// @tparam surface_container the surfaces to be drawn
/// @tparam instance_type the reference instance type
///
/// @param gctx The Geometry context of this operation
/// @param surfaces The surfaces to be converted
/// @param cOptions the covnersion options
/// @param sgi [in,out] the proto indexed grid to be whown
/// @param delegate the delegate to be translated
/// @param refInstance the reference input type from the reference Axes
template <typename surface_container, typename instance_type>
void convertImpl(const GeometryContext& gctx, const surface_container& surfaces,
                 const Options& cOptions, ProtoIndexedSurfaceGrid& sgi,
                 const Experimental::SurfaceCandidatesUpdator& delegate,
                 [[maybe_unused]] const instance_type& refInstance) {
  using GridType =
      typename instance_type::template grid_type<std::vector<std::size_t>>;
  // Defining a Delegate type
  using DelegateType = Experimental::IndexedSurfacesImpl<GridType>;
  // Get the instance
  const auto* instance = delegate.instance();
  auto castedDelegate = dynamic_cast<const DelegateType*>(instance);
  if (castedDelegate != nullptr) {
    auto [pSurfaces, pGrid, pIndices] =
        convert(gctx, surfaces, *castedDelegate, cOptions);
    std::get<0u>(sgi) = pSurfaces;
    std::get<1u>(sgi) = pGrid;
    std::get<2u>(sgi) = pIndices;
  }
}

/// @brief Unrolling function for catching the right instance
///
/// @note parameters are as of the `convertImpl` method
template <typename surface_container, typename tuple_type, std::size_t... I>
void convert(const GeometryContext& gctx, const surface_container& surfaces,
             const Options& cOptions, ProtoIndexedSurfaceGrid& sgi,
             const Experimental::SurfaceCandidatesUpdator& delegate,
             const tuple_type& axesTuple, std::index_sequence<I...>) {
  (convertImpl(gctx, surfaces, cOptions, sgi, delegate, std::get<I>(axesTuple)),
   ...);
}

/// Convert a surface array into needed constituents
///
/// @param gtcx is the geometry context of the conversion call
/// @param surfaces the container of surfaces
/// @param indexGrid the indexGrid delegate
/// @param cOptions the conversion options
///
/// @return a collection of proto surface object and a grid, and associations
template <typename surface_container>
ProtoIndexedSurfaceGrid convertDelegate(
    const GeometryContext& gctx, const surface_container& surfaces,
    const Experimental::SurfaceCandidatesUpdator& delegate,
    const Options& cOptions) {
  // Prep work what is to be filled
  std::vector<ProtoSurface> pSurfaces;
  ProtoGrid pGrid;
  std::vector<std::vector<std::size_t>> indices;
  ProtoIndexedSurfaceGrid sgi = {pSurfaces, pGrid, indices};
  // Convert if dynamic cast happens to work
  convert(
      gctx, surfaces, cOptions, sgi, delegate, Experimental::s_possibleAxes,
      std::make_index_sequence<
          std::tuple_size<decltype(Experimental::s_possibleAxes)>::value>());
  // Return the newly filled ones
  return sgi;
}

/// Convert a surface array into needed constituents
///
/// @param gtcx is the geometry context of the conversion call
/// @param surfaces the container of surfaces
/// @param indexGrid the indexGrid delegate
/// @param cOptions the conversion options
///
/// @return a collection of proto surface object and a grid, and associations
template <typename surface_container, typename index_grid>
ProtoIndexedSurfaceGrid convert(const GeometryContext& gctx,
                                const surface_container& surfaces,
                                const index_grid& indexGrid,
                                const Options& cOptions) {
  // The return surfaces
  std::vector<ProtoSurface> pSurfaces;
  pSurfaces.reserve(surfaces.size());

  for (auto [is, s] : enumerate(surfaces)) {
    // Create the surface converter options
    SurfaceConverter::Options sOptions;
    Style surfaceStyle;
    auto sfIter = cOptions.surfaceStyles.find(s->geometryId());
    if (sfIter != cOptions.surfaceStyles.end()) {
      surfaceStyle = *sfIter;
    }
    sOptions.style = surfaceStyle;
    pSurfaces.push_back(SurfaceConverter::convert(gctx, *s, sOptions));
  }
  // Create the grid
  ProtoGrid pGrid = GridConverter::convert(indexGrid.grid, indexGrid.casts,
                                           cOptions.gridOptions);

  auto axes = indexGrid.grid.axes();

  // Specify the highlight indices
  std::vector<std::vector<std::size_t>> highlightIndices;

  // 1D connections
  if constexpr (index_grid::grid_type::DIM == 1u) {
    for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      typename index_grid::grid_type::index_t lbin;
      lbin[0u] = ib0;
      highlightIndices.push_back(indexGrid.grid.atLocalBins(lbin));
      // Register the bin naming
      pGrid._bin_ids.push_back(std::string("- bin : [") + std::to_string(ib0) +
                               std::string("]"));
    }
  }
  // 2D connections
  if constexpr (index_grid::grid_type::DIM == 2u) {
    for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      for (unsigned int ib1 = 1u; ib1 <= axes[1u]->getNBins(); ++ib1) {
        typename index_grid::grid_type::index_t lbin;
        lbin[0u] = ib0;
        lbin[1u] = ib1;
        highlightIndices.push_back(indexGrid.grid.atLocalBins(lbin));
        // Register the bin naming
        pGrid._bin_ids.push_back(std::string("- bin : [") +
                                 std::to_string(ib0) + std::string(", ") +
                                 std::to_string(ib1) + std::string("]"));
      }
    }
  }
  return std::tie(pSurfaces, pGrid, highlightIndices);
}

}  // namespace IndexedSurfacesConverter

namespace View {

/// Convert into an acts::svg::object with an XY view
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
    if (pGrid._type == actsvg::proto::grid::e_z_phi) {
      sObs.push_back(Acts::Svg::View::zrphi(s, s._name));
    } else {
      sObs.push_back(Acts::Svg::View::xy(s, s._name));
    }
  }

  // The grid
  auto gOb =
      actsvg::display::grid(identification + std::string("_grid"), pGrid);

  // The connectors
  actsvg::connectors::connect_action(gOb._sub_objects, sObs, pIndices);

  // Add them all
  xyIndexedGrid.add_objects(sObs);

  auto xmax = xyIndexedGrid._x_range[1u];
  // The assoication info boxes
  for (auto [ig, gTile] : enumerate(gOb._sub_objects)) {
    // Target surface text
    std::vector<std::string> binText;
    binText.push_back("Source:");
    binText.push_back(pGrid._bin_ids[ig]);
    binText.push_back("Target:");
    for (auto& sis : pIndices[ig]) {
      binText.push_back(std::string("- object: ") + std::to_string(sis));
    }
    // Make the connected text
    auto cText = actsvg::draw::connected_text(
        "bla", {static_cast<actsvg::scalar>(1.1 * xmax), 0}, binText,
        actsvg::style::font{}, actsvg::style::transform{}, gTile);
    xyIndexedGrid.add_object(cText);
  }

  xyIndexedGrid.add_object(gOb);

  // return the grid
  return xyIndexedGrid;
}

}  // namespace View

}  // namespace Svg

}  // namespace Acts