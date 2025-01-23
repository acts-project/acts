// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/detail/IndexedGridFiller.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <algorithm>
#include <array>
#include <memory>

namespace Acts::detail::IndexedSurfacesGenerator {

/// Factory method to create a 1D indexed surface grid
///
/// @param gctx the geometry context
/// @param surfaces the surfaces to be indexed
/// @param rGenerator the reference generator
/// @param pAxis the proto axis
/// @param assignToAll the indices assigned to all bins
/// @param transform the transform into the local binning schema
/// @return an internal navigation delegate
template <typename surface_container, typename reference_generator>
Experimental::InternalNavigationDelegate createInternalNavigation(
    const GeometryContext& gctx, const surface_container& surfaces,
    const reference_generator& rGenerator, const ProtoAxis& pAxis,
    const std::vector<std::size_t> assignToAll = {},
    const Transform3 transform = Transform3::Identity()) {
  // Let the axis create the grid
  return pAxis.getAxis().visit([&]<typename AxisTypeA>(const AxisTypeA& axis)
                                   -> Experimental::InternalNavigationDelegate {
    Experimental::InternalNavigationDelegate nStateUpdater;
    Grid<std::vector<std::size_t>, AxisTypeA> grid(axis);
    std::array<AxisDirection, 1u> axisDirs = {pAxis.getAxisDirection()};

    std::vector<std::size_t> fillExpansion = {pAxis.getFillExpansion()};
    Experimental::detail::IndexedGridFiller filler{fillExpansion};

    // filler.fill(gctx, indexedSurfaces, surfaces, rGenerator, assignToAll);
    // indexed_updator<GridType> indexedSurfaces(std::move(grid), axisDirs,
    //                                           transform);

    return nStateUpdater;
  });
}

/// Factory method to create a 2D indexed surface grid
///
/// @param gctx the geometry context
/// @param surfaces the surfaces to be indexed
/// @param rGenerator the reference generator
/// @param pAxisA the first proto axis
/// @param pAxisB the second proto axis
/// @param assignToAll the indices assigned to all bins
/// @param transform the transform into the local binning schema
///
/// @return an internal navigation delegate
template <typename surface_container, typename reference_generator>
Experimental::InternalNavigationDelegate createInternalNavigation(
    const GeometryContext& gctx, const surface_container& surfaces,
    const reference_generator& rGenerator, const ProtoAxis& pAxisA,
    const ProtoAxis& pAxisB, const std::vector<std::size_t> assignToAll = {},
    const Transform3 transform = Transform3::Identity()) {
  // Let the axes create the grid
  return pAxisA.getAxis().visit(
      [&]<typename AxisTypeA>(
          const AxisTypeA& axisA) -> Experimental::InternalNavigationDelegate {
        return pAxisB.getAxis().visit(
            [&]<typename AxisTypeB>(const AxisTypeB& axisB)
                -> Experimental::InternalNavigationDelegate {
              Experimental::InternalNavigationDelegate nStateUpdater;
              Grid<std::vector<std::size_t>, AxisTypeA, AxisTypeB> grid(axisA,
                                                                        axisB);
              std::array<AxisDirection, 2u> axisDirs = {
                  pAxisA.getAxisDirection(), pAxisB.getAxisDirection()};

              std::vector<std::size_t> fillExpansion = {
                  pAxisA.getFillExpansion(), pAxisB.getFillExpansion()};
              Experimental::detail::IndexedGridFiller filler{fillExpansion};

              // indexed_updator<GridType> indexedSurfaces(std::move(grid),
              // axisDirs,
              //                                           transform);

              return nStateUpdater;
            });
      });
}

}  // namespace Acts::detail::IndexedSurfacesGenerator

/**
  template <typename surface_container>
  InternalNavigationDelegate createInternalNavigation(){

  }

  /// The surfaces to be indexed
  /// (including surfaces that are assigned to all bins)
  surface_container surfaces = {};
  // Indices of surfaces that are to be assigned to all bins
  std::vector<std::size_t> assignToAll = {};
  /// The binning for the indexing
  std::vector<AxisDirection> bValues = {};
  // Bin expansion
  std::vector<std::size_t> binExpansion = {};
  /// The transform into the local binning schema
  Transform3 transform = Transform3::Identity();
  /// Screen output logger
  std::unique_ptr<const Logger> oLogger =
      getDefaultLogger("IndexedSurfacesGenerator", Logging::INFO);

  /// Create the Surface candidate updator
  ///
  /// @tparam axis_generator does generate the axis of the grid
  /// @tparam reference_generator does generate the reference query points
  ///
  /// @param gctx the geometry context
  /// @param aGenerator the axis generator
  /// @param rGenerator the reference generataor
  ///
  /// @return an InternalNavigationDelegate
  template <typename axis_generator, typename reference_generator>
  InternalNavigationDelegate operator()(
      const GeometryContext& gctx, const axis_generator& aGenerator,
      const reference_generator& rGenerator) const {
    ACTS_DEBUG("Indexing " << surfaces.size() << " surface, "
                           << assignToAll.size() << " of which into all bins.");
    // Create the grid with the provided axis generator
    using GridType =
        typename axis_generator::template grid_type<std::vector<std::size_t>>;
    GridType grid(std::move(aGenerator()));

    std::array<AxisDirection, decltype(grid)::DIM> bvArray = {};
    for (auto [ibv, bv] : enumerate(bValues)) {
      bvArray[ibv] = bv;
    }

    indexed_updator<GridType> indexedSurfaces(std::move(grid), bvArray,
                                              transform);
    // Fill the bin indices
    IndexedGridFiller filler{binExpansion};
    filler.oLogger = oLogger->cloneWithSuffix("_filler");
    filler.fill(gctx, indexedSurfaces, surfaces, rGenerator, assignToAll);

    // The portal delegate
    AllPortalsNavigation allPortals;

    // The chained delegate: indexed surfaces and all portals
    using DelegateType =
        IndexedSurfacesAllPortalsNavigation<decltype(grid), indexed_updator>;
    auto indexedSurfacesAllPortals = std::make_unique<const DelegateType>(
        std::tie(allPortals, indexedSurfaces));

    // Create the delegate and connect it
    InternalNavigationDelegate nStateUpdater;
    nStateUpdater.connect<&DelegateType::update>(
        std::move(indexedSurfacesAllPortals));
    return nStateUpdater;
  }

  /// Access to the logger
  ///
  /// @return a const reference to the logger
  const Logger& logger() const { return *oLogger; }
};

}  // namespace Acts::Experimental::detail
*/