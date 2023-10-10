// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Detector/detail/IndexedSurfacesGenerator.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Plugins/ActSVG/GridSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/detail/Grid.hpp"
#include <actsvg/core.hpp>
#include <actsvg/meta.hpp>

#include <tuple>
#include <vector>

namespace Acts {

namespace Experimental {

using namespace detail::GridAxisGenerators;

// Generate the possible axes in this case
static auto s_possibleAxes =
    std::tuple<EqBound, EqOpen, EqClosed,
               // All 1D Var  options
               VarBound, VarOpen, VarClosed,
               // All 2D EqEq options
               EqBoundEqBound, EqBoundEqOpen, EqBoundEqClosed, EqOpenEqBound,
               EqOpenEqOpen, EqOpenEqClosed, EqClosedEqBound, EqClosedEqOpen,
               EqClosedEqClosed,
               // All 2D EqVar options
               EqBoundVarBound, EqBoundVarOpen, EqBoundVarClosed,
               EqOpenVarBound, EqOpenVarOpen, EqOpenVarClosed, EqClosedVarBound,
               EqClosedVarOpen, EqClosedVarClosed,
               // All 2D VarEq options
               VarBoundEqBound, VarBoundEqOpen, VarBoundEqClosed,
               VarOpenEqBound, VarOpenEqOpen, VarOpenEqClosed, VarClosedEqBound,
               VarClosedEqOpen, VarClosedEqClosed,
               // All 2D VarEq options
               VarBoundVarBound, VarBoundVarOpen, VarBoundVarClosed,
               VarOpenVarBound, VarOpenVarOpen, VarOpenVarClosed,
               VarClosedVarBound, VarClosedVarOpen, VarClosedVarClosed>{};

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
};

/// Convert a surface array into needed constituents
///
/// @note actual conversion implementation, bottom of unrolling loop
///
/// @param gtcx is the geometry context of the conversion call
/// @param surfaces the container of surfaces
/// @param indexGrid the indexGrid delegate
/// @param cOptions the conversion options
///
/// @return a collection of proto surface object and a grid, and associations
template <typename surface_container, typename index_grid>
ProtoIndexedSurfaceGrid convertImpl(const GeometryContext& gctx,
                                    const surface_container& surfaces,
                                    const index_grid& indexGrid,
                                    const Options& cOptions) {
  // The return surfaces
  std::vector<ProtoSurface> pSurfaces;
  pSurfaces.reserve(surfaces.size());

  // Make a local copy of the grid Options, in case we have to estimate
  // an additional bound
  GridConverter::Options gridOptions = cOptions.gridOptions;
  Extent constrain;

  // Estimate the radial extension
  // - for 1D phi
  // - for 2D z-phi or phi-z
  bool estimateR =
      (index_grid::grid_type::DIM == 1 and indexGrid.casts[0u] == binPhi) or
      (index_grid::grid_type::DIM == 2 and
       (indexGrid.casts[0u] == binPhi or indexGrid.casts[1u] == binPhi));

  for (auto [is, s] : enumerate(surfaces)) {
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
        pSurface._radii[0u] = sExtent.medium(binR);
      }
      constrain.extend(sExtent, {binR});
    }
    pSurfaces.push_back(pSurface);
  }

  // Adjust the grid options
  if constexpr (index_grid::grid_type::DIM == 1u) {
    if (indexGrid.casts[0u] == binPhi) {
      auto estRangeR = constrain.range(binR);
      std::array<ActsScalar, 2u> rRange = {estRangeR.min(), estRangeR.max()};
      gridOptions.optionalBound = {rRange, binR};
    }
  }

  // Create the grid
  ProtoGrid pGrid =
      GridConverter::convert(indexGrid.grid, indexGrid.casts, gridOptions);

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
        if (estimateR) {
          pGrid._reference_r = constrain.medium(binR);
        }
      }
    }
  }
  return std::tie(pSurfaces, pGrid, highlightIndices);
}

/// @brief Convert the single delegate if it is of the type of the reference
///
/// @note It will do nothing if the type does not match
///
/// @tparam surface_container the surfaces to be drawn
/// @tparam instance_type the reference instance type
///
/// @param gctx The Geometry context of this operation
/// @param surfaces The surfaces to be converted
/// @param cOptions the conversion options
/// @param sgi [in,out] the proto indexed grid to be converted
/// @param delegate the delegate to be translated
/// @param refInstance the reference input type from the reference Axes
template <typename surface_container, typename instance_type>
void convert(const GeometryContext& gctx, const surface_container& surfaces,
             const Options& cOptions, ProtoIndexedSurfaceGrid& sgi,
             const Experimental::SurfaceCandidatesUpdator& delegate,
             [[maybe_unused]] const instance_type& refInstance) {
  using GridType =
      typename instance_type::template grid_type<std::vector<std::size_t>>;
  // Defining a Delegate type
  using DelegateType = Experimental::IndexedSurfacesAllPortalsImpl<
      GridType, Experimental::IndexedSurfacesImpl>;
  using SubDelegateType = Experimental::IndexedSurfacesImpl<GridType>;

  // Get the instance
  const auto* instance = delegate.instance();
  auto castedDelegate = dynamic_cast<const DelegateType*>(instance);
  if (castedDelegate != nullptr) {
    // Get the surface updator
    auto indexedSurfaces = std::get<SubDelegateType>(castedDelegate->updators);
    auto [pSurfaces, pGrid, pIndices] =
        convertImpl(gctx, surfaces, indexedSurfaces, cOptions);
    std::get<0u>(sgi) = pSurfaces;
    std::get<1u>(sgi) = pGrid;
    std::get<2u>(sgi) = pIndices;
  }
}

/// @brief Unrolling function for catching the right instance
///
/// @note parameters are as of the `convertImpl` method
template <typename surface_container, typename tuple_type, std::size_t... I>
void unrollConvert(const GeometryContext& gctx,
                   const surface_container& surfaces, const Options& cOptions,
                   ProtoIndexedSurfaceGrid& sgi,
                   const Experimental::SurfaceCandidatesUpdator& delegate,
                   const tuple_type& axesTuple, std::index_sequence<I...>) {
  (convert(gctx, surfaces, cOptions, sgi, delegate, std::get<I>(axesTuple)),
   ...);
}

/// Convert a surface array into needed constituents
///
/// @param gtcx is the geometry context of the conversion call
/// @param surfaces the container of surfaces
/// @param indexGrid the indexGrid delegate
/// @param cOptions the conversion options
///
/// @note this is the entry point of the conversion, i.e. top of the
/// unrolling loop
///
/// @return a collection of proto surface object and a grid, and associations
template <typename surface_container>
ProtoIndexedSurfaceGrid convert(
    const GeometryContext& gctx, const surface_container& surfaces,
    const Experimental::SurfaceCandidatesUpdator& delegate,
    const Options& cOptions) {
  // Prep work what is to be filled
  std::vector<ProtoSurface> pSurfaces;
  ProtoGrid pGrid;
  std::vector<std::vector<std::size_t>> indices;
  ProtoIndexedSurfaceGrid sgi = {pSurfaces, pGrid, indices};
  // Convert if dynamic cast happens to work
  unrollConvert(
      gctx, surfaces, cOptions, sgi, delegate, Experimental::s_possibleAxes,
      std::make_index_sequence<
          std::tuple_size<decltype(Experimental::s_possibleAxes)>::value>());
  // Return the newly filled ones
  return sgi;
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
  // The association info boxes
  for (auto [ig, gTile] : enumerate(gOb._sub_objects)) {
    // Target surface text
    std::vector<std::string> binText;
    binText.push_back("Source:");
    binText.push_back(pGrid._bin_ids[ig]);
    binText.push_back("Target:");
    for (const auto& sis : pIndices[ig]) {
      binText.push_back(std::string("- object: ") + std::to_string(sis));
    }
    // Make the connected text
    std::string cTextId =
        identification + std::string("_ct_") + std::to_string(ig);
    auto cText = actsvg::draw::connected_text(
        "bla", {static_cast<actsvg::scalar>(1.1 * xmax), 0}, binText,
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
}  // namespace Svg
}  // namespace Acts
