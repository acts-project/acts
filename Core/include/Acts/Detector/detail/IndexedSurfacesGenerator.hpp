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
template <template <typename> class indexed_updator, typename surface_container,
          typename reference_generator>
Experimental::InternalNavigationDelegate createInternalNavigation(
    const GeometryContext& gctx, const surface_container& surfaces,
    const reference_generator& rGenerator, const ProtoAxis& pAxis,
    const std::vector<std::size_t> assignToAll = {},
    const Transform3 transform = Transform3::Identity()) {
  // Let the axis create the grid
  return pAxis.getAxis().visit([&]<typename AxisTypeA>(const AxisTypeA& axis) {
    Grid<std::vector<std::size_t>, AxisTypeA> grid(axis);

    // Prepare the indexed updator
    std::array<AxisDirection, 1u> axisDirs = {pAxis.getAxisDirection()};
    indexed_updator<decltype(grid)> indexedSurfaces(std::move(grid), axisDirs,
                                                    transform);

    // Prepare the filling
    Experimental::InternalNavigationDelegate nStateUpdater;

    std::vector<std::size_t> fillExpansion = {pAxis.getFillExpansion()};
    Experimental::detail::IndexedGridFiller filler{fillExpansion};
    filler.fill(gctx, indexedSurfaces, surfaces, rGenerator, assignToAll);

    // The portal delegate
    Experimental::AllPortalsNavigation allPortals;

    // The chained delegate: indexed surfaces and all portals
    using DelegateType =
        Experimental::IndexedSurfacesAllPortalsNavigation<decltype(grid),
                                                          indexed_updator>;
    auto indexedSurfacesAllPortals = std::make_unique<const DelegateType>(
        std::tie(allPortals, indexedSurfaces));

    // Create the delegate and connect it
    nStateUpdater.connect<&DelegateType::update>(
        std::move(indexedSurfacesAllPortals));

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
template <template <typename> class indexed_updator, typename surface_container,
          typename reference_generator>
Experimental::InternalNavigationDelegate createInternalNavigation(
    const GeometryContext& gctx, const surface_container& surfaces,
    const reference_generator& rGenerator, const ProtoAxis& pAxisA,
    const ProtoAxis& pAxisB, const std::vector<std::size_t> assignToAll = {},
    const Transform3 transform = Transform3::Identity()) {
  // Let the axes create the grid
  return pAxisA.getAxis().visit([&]<typename AxisTypeA>(
                                    const AxisTypeA& axisA) {
    return pAxisB.getAxis().visit([&]<typename AxisTypeB>(
                                      const AxisTypeB& axisB) {
      Grid<std::vector<std::size_t>, AxisTypeA, AxisTypeB> grid(axisA, axisB);
      Experimental::InternalNavigationDelegate nStateUpdater;

      // Prepare the indexed updator
      std::array<AxisDirection, 2u> axisDirs = {pAxisA.getAxisDirection(),
                                                pAxisB.getAxisDirection()};
      indexed_updator<decltype(grid)> indexedSurfaces(std::move(grid), axisDirs,
                                                      transform);

      std::vector<std::size_t> fillExpansion = {pAxisA.getFillExpansion(),
                                                pAxisB.getFillExpansion()};

      Experimental::detail::IndexedGridFiller filler{fillExpansion};
      filler.fill(gctx, indexedSurfaces, surfaces, rGenerator, assignToAll);

      // The portal delegate
      Experimental::AllPortalsNavigation allPortals;

      // The chained delegate: indexed surfaces and all portals
      using DelegateType =
          Experimental::IndexedSurfacesAllPortalsNavigation<decltype(grid),
                                                            indexed_updator>;
      auto indexedSurfacesAllPortals = std::make_unique<const DelegateType>(
          std::tie(allPortals, indexedSurfaces));

      // Create the delegate and connect it
      nStateUpdater.connect<&DelegateType::update>(
          std::move(indexedSurfacesAllPortals));

      return nStateUpdater;
    });
  });
}

}  // namespace Acts::detail::IndexedSurfacesGenerator
