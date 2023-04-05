// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Utilities/Enumerate.hpp"

namespace Acts {
namespace Experimental {

/// @brief  A templaed indexed grid generator.
///
/// This Generator creates a SurfaceCandidatesUpdator delegate
/// which can then be used in the DetectorVolume class for updating
/// given surface candidates.
///
/// @tparam objects_container the objects container
template <typename surface_container>
struct IndexedSurfacesGenerator {
  /// The surfaces to be indexed
  /// (including surfaces that are assigned to all bins)
  surface_container surfaces = {};
  // Indices of surfaces that are to be assigned to all bins
  std::vector<std::size_t> assignToAll = {};
  /// The binning
  std::vector<BinningValue> bValues = {};

  // Bin expansion
  std::vector<std::size_t> binExpansion = {};

  /// The transform into the local binning schema
  Transform3 transform = Transform3::Identity();

  // An ouput logging level
  Logging::Level logLevel = Logging::INFO;

  /// Create the Surface candidate updator
  ///
  /// @tparam axis_generator does generate the axis of the grid
  /// @tparam reference_generator does generate the reference query points
  ///
  /// @param gctx the geometry context
  /// @param aGenerator the axis generator
  /// @param rGenerator the reference generataor
  ///
  /// @return a SurfaceCandidateUpdator delegate
  template <typename axis_generator, typename reference_generator>
  SurfaceCandidatesUpdator operator()(
      const GeometryContext& gctx, const axis_generator& aGenerator,
      const reference_generator& rGenerator) const {
    ACTS_LOCAL_LOGGER(getDefaultLogger("IndexedSurfacesGenerator", logLevel));
    ACTS_DEBUG("Indexing " << surfaces.size() << " surface, "
                           << assignToAll.size() << " of which into all bins.");
    // Create the grid with the provided axis generator
    typename axis_generator::template grid_type<std::vector<std::size_t>> grid(
        std::move(aGenerator()));

    std::array<BinningValue, decltype(grid)::DIM> bvArray = {};
    for (auto [ibv, bv] : enumerate(bValues)) {
      bvArray[ibv] = bv;
    }

    // Create the indexed surface grid:
    // non-const as we need to fill the grid
    auto indexedSurfaces =
        std::make_unique<IndexedSurfacesImpl<decltype(grid)>>(
            std::move(grid), bvArray, transform);

    // Fill the bin indicies
    IndexedGridFiller filler{binExpansion, logLevel};
    filler.fill(gctx, *indexedSurfaces, surfaces, rGenerator, assignToAll);

    // Now turn it const to act as a navigation delegate
    std::unique_ptr<const IndexedSurfacesImpl<decltype(grid)>> cis =
        std::move(indexedSurfaces);

    // Now create the delegate
    SurfaceCandidatesUpdator nStateUpdator;
    nStateUpdator.connect<&IndexedSurfacesImpl<decltype(grid)>::update>(
        std::move(cis));
    return nStateUpdator;
  }
};

}  // namespace Experimental
}  // namespace Acts
