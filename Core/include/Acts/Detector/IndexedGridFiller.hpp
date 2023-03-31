// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <set>
#include <string>
#include <vector>

namespace {

/// @brief Helper method to loop correctly in local bins given open, bound, closed
///
/// @param minMaxBins estimated
/// @param expand the parameter to expand the view ( extra window)
/// @param nBins the maximum number of bins in this
/// @param type the boundary type of the axis
///
/// @note for closed binning a span over half the bins flips direction
///
/// @return a vector of bins to be filled
std::vector<std::size_t> binSequence(std::array<std::size_t, 2u> minMaxBins,
                                     std::size_t expand, std::size_t nBins,
                                     Acts::detail::AxisBoundaryType type) {
  // Return vector for iterations
  std::vector<std::size_t> rBins;
  /// Helper method to fill a range
  auto fill_linear = [&](std::size_t lmin, std::size_t lmax) -> void {
    for (std::size_t b = lmin; b <= lmax; ++b) {
      rBins.push_back(b);
    }
  };
  std::size_t bmin = minMaxBins[0u];
  std::size_t bmax = minMaxBins[1u];

  // Open/Bound cases
  if (type != Acts::detail::AxisBoundaryType::Closed) {
    rBins.reserve(bmax - bmin + 1u + 2 * expand);
    // handle bmin:/max expand it down (for bound, don't fill underflow)
    if (type == Acts::detail::AxisBoundaryType::Bound) {
      bmin = (int(bmin) - int(expand) > 0) ? bmin - expand : 1u;
      bmax = (bmax + expand <= nBins) ? bmax + expand : nBins;
    } else if (type == Acts::detail::AxisBoundaryType::Open) {
      bmin = (int(bmin) - int(expand) >= 0u) ? bmin - expand : 0u;
      bmax = (bmax + expand <= nBins + 1u) ? bmax + expand : nBins + 1u;
    }
    fill_linear(bmin, bmax);
  } else {
    std::size_t span = bmax - bmin + 1u + 2 * expand;
    // Safe with respect to the closure point, treat as bound
    if (2 * span < nBins and (bmax + expand <= nBins) and
        (int(bmin) - int(expand) > 0)) {
      return binSequence({bmin, bmax}, expand, nBins,
                         Acts::detail::AxisBoundaryType::Bound);
    } else if (2 * span < nBins) {
      bmin = int(bmin) - int(expand) > 0 ? bmin - expand : 1u;
      bmax = bmax + expand <= nBins ? bmax + expand : nBins;
      fill_linear(bmin, bmax);
      // deal with expansions over the phi boundary
      if (bmax + expand > nBins) {
        std::size_t overstep = (bmax + expand - nBins);
        fill_linear(1u, overstep);
      }
      if (int(bmin) - int(expand) < 1) {
        std::size_t understep = abs(int(bmin) - int(expand));
        fill_linear(nBins - understep, nBins);
      }
      std::sort(rBins.begin(), rBins.end());
    } else {
      // Jump over the phi boundary
      fill_linear(bmax - expand, nBins);
      fill_linear(1, bmin + expand);
      std::sort(rBins.begin(), rBins.end());
    }
  }
  return rBins;
}

/// @brief Helper method to fill intermediate bins in a given set
///
/// @tparam grid_type the type of the grid that determines locall binning
///
/// @param grid the grid used for this
/// @param queries the grid positions for the test
/// @param expand the expansion (symmetric) around the bins
template <typename grid_type>
std::set<typename grid_type::index_t> localIndices(
    const grid_type& grid,
    const std::vector<typename grid_type::point_t>& queries,
    std::size_t expand = 0u) {
  // Return indices
  std::set<typename grid_type::index_t> lIndices;

  if (queries.empty()) {
    throw std::runtime_error("IndexedSurfaceGridFiller: no query point given.");
  }

  /// These are the axis bounds type parameters - for correct bin sequences
  std::array<Acts::detail::AxisBoundaryType, grid_type::DIM> axisTypes{};
  std::array<std::size_t, grid_type::DIM> axisBins{};
  // Fill the axis types
  for (auto [ia, a] : enumerate(grid.axes())) {
    axisTypes[ia] = a->getBoundaryType();
    axisBins[ia] = a->getNBins();
  }

  // Initialize the bin ranges
  std::array<std::array<std::size_t, 2u>, grid_type::DIM> binRanges;
  for (auto& br : binRanges) {
    br[0u] = std::numeric_limits<std::size_t>::max();
    br[1u] = 0u;
  }
  // Bin range bounding box - estimated from the query points
  for (const auto& q : queries) {
    auto qbin = grid.localBinsFromPosition(q);
    for (std::size_t ib = 0; ib < grid_type::DIM; ++ib) {
      auto iqb = qbin[ib];
      binRanges[ib][0u] = iqb < binRanges[ib][0u] ? iqb : binRanges[ib][0u];
      binRanges[ib][1u] = iqb > binRanges[ib][1u] ? iqb : binRanges[ib][1u];
    }
  }
  // Fill the bins - 1D case
  if constexpr (grid_type::DIM == 1u) {
    auto localBins0 =
        binSequence(binRanges[0u], expand, axisBins[0u], axisTypes[0u]);
    for (auto l0 : localBins0) {
      typename grid_type::index_t b;
      b[0u] = l0;
      lIndices.insert(b);
    }
  }
  // Fill the bins - 2D case
  if constexpr (grid_type::DIM == 2u) {
    auto localBins0 =
        binSequence(binRanges[0u], expand, axisBins[0u], axisTypes[0u]);
    auto localBins1 =
        binSequence(binRanges[1u], expand, axisBins[1u], axisTypes[1u]);
    for (auto l0 : localBins0) {
      for (auto l1 : localBins1) {
        typename grid_type::index_t b;
        b[0u] = l0;
        b[1u] = l1;
        lIndices.insert(b);
      }
    }
  }
  return lIndices;
}

/// @brief Helper method to screen output the local bins
///
/// @tparam local_bin the type of the local bins
///
/// @param lbins the local bins
///
/// @return a string containing the local bins orderd in a set
template <typename local_bin>
std::string outputIndices(const std::set<local_bin>& lbins) {
  std::string rString;
  for (auto [ilb, lb] : Acts::enumerate(lbins)) {
    if (ilb == 0) {
      rString = "bins: [";
    } else {
      rString += ", [";
    }
    for (auto [ib, b] : Acts::enumerate(lb)) {
      if (ib != 0u) {
        rString += ", ";
      }
      rString += std::to_string(b);
    }
    rString += "]";
  }
  return rString;
}

}  // namespace

namespace Acts {
namespace Experimental {

/// A struct to access the center position
///
/// This generator will provide only one filling point and hence
/// only a single bin in the indexed grid.
struct CenterReferenceGenerator {
  /// Helper to access the Center point of for filling the grid
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference point is to be accessed
  ///
  /// @return a vector of referene points for filling
  const std::vector<Vector3> references(const GeometryContext& gctx,
                                        const Surface& surface) const {
    return {surface.center(gctx)};
  }
};

/// A struct to access reference postions based on bin values
///
/// This generator will provide only one filling point and hence
/// only a single bin in the indexed grid.
struct BinningValueReferenceGenerator {
  /// The binning value
  BinningValue bValue;

  /// Helper to access a reference postion based on binning value
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference point is to be accessed
  ///
  /// @return a vector of referene points for filling
  const std::vector<Vector3> references(const GeometryContext& gctx,
                                        const Surface& surface) const {
    return {surface.binningPosition(gctx, bValue)};
  }
};

/// A struct to access generated vertices from surface polyhedrons
/// These vertices are then used to find the bin boundary box for the
/// indexed grid.
///
/// The grid filling then completes the empty bins in between and
/// expands if necessary.
struct PolyhedronReferenceGenerator {
  /// Also use the barycenter
  bool addBarycenter = true;

  /// The number of segments to approximate (1 means extrema points only)
  unsigned int nSegments = 1;

  /// Helper to access the Center point of for filling the grid
  ///
  /// @param gctx the geometry context of this operation
  /// @param surface the surface for which the reference point is to be accessed
  ///
  /// @return a vector of referene points for filling
  const std::vector<Vector3> references(const GeometryContext& gctx,
                                        const Surface& surface) const {
    // Create the return  vector
    std::vector<Vector3> rPositions;
    auto pHedron = surface.polyhedronRepresentation(gctx, nSegments);
    rPositions.insert(rPositions.end(), pHedron.vertices.begin(),
                      pHedron.vertices.end());
    // Add the barycenter if configured
    if (addBarycenter) {
      Vector3 bc(0., 0., 0.);
      std::for_each(rPositions.begin(), rPositions.end(),
                    [&](const auto& p) { bc += p; });
      bc *= 1. / rPositions.size();
      rPositions.push_back(bc);
    }
    return rPositions;
  }
};

/// A helper class that fills surfaces into predefined grids
struct IndexedGridFiller {
  // Bin expansion where needed
  std::size_t binExpansion = 0u;

  // An ouput logging level
  Logging::Level logLevel = Logging::INFO;

  /// @brief This method takes a collection of objects and fills them
  /// into an index grid
  ///
  /// @tparam index_grid the type of the index grid
  /// @tparam indexed_objects the type of the object container
  /// @tparam reference_generator the generator for reference points to be filled
  ///
  /// @param gctx the geometry context of the operation
  /// @param iGrid [in,out] the index grid object to be filled
  /// @param iObjects the object container to be indexed
  /// @param rGenerator the reference point generator for position queries
  ///
  /// @note as this is a Detector module, the objects within the indexed_objects container
  /// are assumed to have pointer semantics
  ///
  template <typename index_grid, typename indexed_objects,
            typename reference_generator>
  void fill(const GeometryContext& gctx, index_grid& iGrid,
            const indexed_objects& iObjects,
            const reference_generator& rGenerator) {
    // Screen output
    ACTS_LOCAL_LOGGER(getDefaultLogger("IndexedGridFiller", logLevel));
    // Loop over the surfaces to be filled
    for (auto [io, o] : enumerate(iObjects)) {
      // Get the reference positions
      auto refs = rGenerator.references(gctx, *o);
      std::vector<typename index_grid::grid_type::point_t> gridQueries;
      gridQueries.reserve(refs.size());
      for (const auto ref : refs) {
        // Cast the transfrom according to the grid binning
        gridQueries.push_back(iGrid.castPosition(ref));
      }
      ACTS_DEBUG(gridQueries.size() << " reference points generated.");
      auto lIndices = localIndices<decltype(iGrid.grid)>(
          iGrid.grid, gridQueries, binExpansion);
      ACTS_DEBUG(lIndices.size() << " indices assigned.");
      if (logLevel <= Logging::VERBOSE) {
        ACTS_VERBOSE("- list of indices: " << outputIndices(lIndices));
      }
      // Now fill the surface indices
      for (const auto& li : lIndices) {
        auto& bContent = iGrid.grid.atLocalBins(li);
        if (std::find(bContent.begin(), bContent.end(), io) == bContent.end()) {
          bContent.push_back(io);
        }
      }
    }
  }
};

}  // namespace Experimental
}  // namespace Acts
