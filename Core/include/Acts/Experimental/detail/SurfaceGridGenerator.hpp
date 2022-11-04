// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/NavigationDelegates.hpp"
#include "Acts/Experimental/detail/NavigationStateUpdators.hpp"
#include "Acts/Experimental/detail/SurfaceGridHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <exception>
#include <optional>
#include <vector>

namespace Acts {

class Surface;

namespace Experimental {
namespace detail {

using AxisEqB = Acts::detail::Axis<Acts::detail::AxisType::Equidistant,
                                   Acts::detail::AxisBoundaryType::Bound>;
using AxisEqC = Acts::detail::Axis<Acts::detail::AxisType::Equidistant,
                                   Acts::detail::AxisBoundaryType::Closed>;

using AxisVarB = Acts::detail::Axis<Acts::detail::AxisType::Variable,
                                    Acts::detail::AxisBoundaryType::Bound>;
using AxisVarC = Acts::detail::Axis<Acts::detail::AxisType::Variable,
                                    Acts::detail::AxisBoundaryType::Closed>;

using IndexContainer = std::vector<std::size_t>;

// The possible grids for surface ordering
using GridEqBEqB = Acts::detail::Grid<IndexContainer, AxisEqB, AxisEqB>;

using GridEqBEqC = Acts::detail::Grid<IndexContainer, AxisEqB, AxisEqC>;

using GridEqBVarB = Acts::detail::Grid<IndexContainer, AxisEqB, AxisVarB>;

using GridEqBVarC = Acts::detail::Grid<IndexContainer, AxisEqB, AxisVarC>;

using GridVarBVarB = Acts::detail::Grid<IndexContainer, AxisVarB, AxisVarB>;

using GridVarBEqB = Acts::detail::Grid<IndexContainer, AxisVarB, AxisEqB>;

using GridVarBEqC = Acts::detail::Grid<IndexContainer, AxisVarB, AxisEqC>;

using GridVarBVarC = Acts::detail::Grid<IndexContainer, AxisVarB, AxisVarC>;

// The possible grid attachers for surface ordering
using AttacherEqBEqB = GridSurfaceAttacher<GridEqBEqB>;
using AttacherEqBEqC = GridSurfaceAttacher<GridEqBEqC>;
using AttacherEqBVarB = GridSurfaceAttacher<GridEqBVarB>;
using AttacherEqBVarC = GridSurfaceAttacher<GridEqBVarC>;
using AttacherVarBVarB = GridSurfaceAttacher<GridVarBVarB>;
using AttacherVarBEqB = GridSurfaceAttacher<GridVarBEqB>;
using AttacherVarBEqC = GridSurfaceAttacher<GridVarBEqC>;
using AttacherVarBVarC = GridSurfaceAttacher<GridVarBVarC>;

auto portalAttacher = AllPortalsAttacher{};

// The possible surface grid Attachers
using UpdatorEqBEqB =
    NavigationStateUpdatorImpl<AttacherEqBEqB, decltype(portalAttacher)>;

using UpdatorEqBEqC =
    NavigationStateUpdatorImpl<AttacherEqBEqC, decltype(portalAttacher)>;

using UpdatorEqBVarB =
    NavigationStateUpdatorImpl<AttacherEqBVarB, decltype(portalAttacher)>;

using UpdatorEqBVarC =
    NavigationStateUpdatorImpl<AttacherEqBVarC, decltype(portalAttacher)>;

using UpdatorVarBVarB =
    NavigationStateUpdatorImpl<AttacherVarBVarB, decltype(portalAttacher)>;

using UpdatorVarBEqB =
    NavigationStateUpdatorImpl<AttacherVarBEqB, decltype(portalAttacher)>;

using UpdatorVarBEqC =
    NavigationStateUpdatorImpl<AttacherVarBEqC, decltype(portalAttacher)>;

using UpdatorVarBVarC =
    NavigationStateUpdatorImpl<AttacherVarBVarC, decltype(portalAttacher)>;

// Possible types
static std::tuple<UpdatorEqBEqB, UpdatorEqBEqC, UpdatorEqBVarB, UpdatorEqBVarC,
                  UpdatorVarBVarB, UpdatorVarBEqB, UpdatorVarBEqC,
                  UpdatorVarBVarC>
    s_gridUpdatorTempaltes = {
        UpdatorEqBEqB(std::make_tuple(
            AttacherEqBEqB(GridEqBEqB(std::make_tuple(AxisEqB(0., 1., 1u),
                                                      AxisEqB(0., 1., 1u))),
                           {binX, binY}),
            portalAttacher)),
        UpdatorEqBEqC(std::make_tuple(
            AttacherEqBEqC(GridEqBEqC(std::make_tuple(AxisEqB(0., 1., 1u),
                                                      AxisEqC(0., 1., 1u))),
                           {binX, binY}),
            portalAttacher)),
        UpdatorEqBVarB(std::make_tuple(
            AttacherEqBVarB(
                GridEqBVarB(std::make_tuple(AxisEqB(0., 1., 1u), AxisVarB({}))),
                {binX, binY}),
            portalAttacher)),
        UpdatorEqBVarC(std::make_tuple(
            AttacherEqBVarC(
                GridEqBVarC(std::make_tuple(AxisEqB(0., 1., 1u), AxisVarC({}))),
                {binX, binY}),
            portalAttacher)),
        UpdatorVarBVarB(std::make_tuple(
            AttacherVarBVarB(
                GridVarBVarB(std::make_tuple(AxisVarB({}), AxisVarB({}))),
                {binX, binY}),
            portalAttacher)),
        UpdatorVarBEqB(std::make_tuple(
            AttacherVarBEqB(
                GridVarBEqB(std::make_tuple(AxisVarB({}), AxisEqB(0., 1., 1u))),
                {binX, binY}),
            portalAttacher)),
        UpdatorVarBEqC(std::make_tuple(
            AttacherVarBEqC(
                GridVarBEqC(std::make_tuple(AxisVarB({}), AxisEqC(0., 1., 1u))),
                {binX, binY}),
            portalAttacher)),
        UpdatorVarBVarC(std::make_tuple(
            AttacherVarBVarC(
                GridVarBVarC(std::make_tuple(AxisVarB({}), AxisVarC({}))),
                {binX, binY}),
            portalAttacher))};

// Helpers to extract grid parameters

// Grid parameters for further digestion
struct AxisParameters {
  const Acts::IAxis* raw = nullptr;
  std::vector<ActsScalar> edges = {};
  Acts::detail::AxisBoundaryType boundaryType =
      Acts::detail::AxisBoundaryType::Closed;
  BinningValue bValue = binX;
};

struct GridParameters {
  // The collected information
  std::vector<AxisParameters> axes;
  std::vector<std::size_t> globalBins;
  std::vector<std::vector<std::size_t>> entries;
};

/// @brief Helper function to extract grid parameter types for drawing/IO
///
/// @tparam T the type of the delegate for checking
/// @param gp the grid parameters vector
/// @param d the delegate in question
/// @param t the reference object for the dynamic_cast
/// @param overflow take also the overflow binds
///
template <typename T>
void extract(GridParameters& gp, const Experimental::IDelegateImpl& d,
             [[maybe_unused]] const T& t, [[maybe_unused]] bool overflow) {
  auto dct = dynamic_cast<const T*>(&d);
  if (dct != nullptr) {
    // Retrieve attacher, grid and casts
    auto& attachers = dct->candidateAttachers;
    auto& gAttacher = std::get<0u>(attachers);
    auto& grid = gAttacher.grid;
    auto& casts = gAttacher.casts;
    const auto axes = grid.axes();


    std::vector<std::vector<ActsScalar>> gridBinCenters = {};
    // Collect the grid boundaries
    for (auto [ia, axis] : enumerate(axes)) {
      // Prepare the axis parameters
      AxisParameters ap;
      ap.raw = axis;
      ap.edges = axis->getBinEdges();
      ap.boundaryType = axis->getBoundaryType();
      ap.bValue = casts[ia];
      gp.axes.push_back(ap);
      // The bin centers for content lookup
      // -> This could go into a helper function, indeed
      std::vector<ActsScalar> axisBinCenters;
      axisBinCenters.reserve(ap.edges.size());
      for (auto [ie, e] : enumerate(ap.edges)) {
        if (ie > 0u) {
          axisBinCenters.push_back(0.5 * (e + ap.edges[ie - 1u]));
        }
      }
      gridBinCenters.push_back(axisBinCenters);
    }

    // There must be a better way to do this generic -------------
    std::size_t nBins = grid.size(false);
    std::vector<std::array<ActsScalar, axes.size()>> gridLookups(
        nBins, std::array<ActsScalar, axes.size()>());

    std::size_t b01 = 0u;
    for (const auto& b0 : gridBinCenters[0u]) {
      for (const auto& b1 : gridBinCenters[1u]) {
        auto& blu = gridLookups[b01++];
        blu[0u] = b0;
        blu[1u] = b1;
      }
    }
    // now look up an dfill
    for (const auto& glue : gridLookups) {
      gp.globalBins.push_back(grid.globalBinFromPosition(glue));
      gp.entries.push_back(grid.atPosition(glue));
    }
    // -----------------------------------------------------------
  }
}

/// Hepler method to extract grid parameters from a delegate
/// @param d the navigation delegate that is being casted
/// @param grid the possible delegate types
/// @param overflow take also the overflow bins
///
/// @return a reolved parameter object
template <typename Grids>
GridParameters extractGridParameters(const Experimental::IDelegateImpl& d,
                                     const Grids& grids,
                                     [[maybe_unused]] bool overflow = false) {
  GridParameters gp;
  std::apply([&](auto&&... g) { ((extract(gp, d, g, overflow)), ...); }, grids);
  return gp;
}

// Avoid doing work twice, i.e. store the polyhedrons
using SurfacPolyhedron = std::tuple<const Surface*, Polyhedron>;

// Poor man's clusterization - to be replaced -------
struct Cluster {
  /// A cluster has a center
  ActsScalar center = 0.;
  /// And remembers its values
  std::vector<ActsScalar> values = {};

  /// Update the cluster with this value
  /// @param v the value to be added
  void update(const ActsScalar& v) {
    values.push_back(v);
    center = 0.;
    std::for_each(values.begin(), values.end(),
                  [&](const auto& v) { center += v; });
    center /= values.size();
  }
};

/// @brief  Poor man's clustering
///
/// @param values
/// @param tolerance
///
/// @return a list of clusters
std::vector<Cluster> clusterize(const std::vector<ActsScalar>& values,
                                const ActsScalar& tolerance) {
  std::vector<Cluster> clusters;
  for (const auto& v : values) {
    bool found = false;
    for (auto& c : clusters) {
      if (std::abs(c.center - v) < tolerance) {
        c.update(v);
        found = true;
        break;
      }
    }
    if (not found) {
      clusters.push_back(Cluster{v, {v}});
    }
  }
  return clusters;
}

// -------------------------

using AxisDefsEq = std::tuple<ActsScalar, ActsScalar, std::size_t,
                              Acts::detail::AxisBoundaryType>;
using AxisDefsVar =
    std::tuple<std::vector<ActsScalar>, Acts::detail::AxisBoundaryType>;
struct AxisDefs {
  std::optional<AxisDefsEq> eqAxisDefs = std::nullopt;
  std::optional<AxisDefsVar> varAxisDefs = std::nullopt;
};

// Get the equidistant parameters
///
/// @param clusters [in,out] the cluster values
/// @param deltaClusters [in,out] the cluster values
/// @param deltaTolerance the toelrance for equidistant
/// @param bValue the binning value
/// @param phiCenter is a center phi value
/// @param forceEquidistant if there are only that amount of values
///
/// @note clusters are extected to be sorted
AxisDefs axisDefinition(std::vector<Cluster>& clusters,
                        ActsScalar diffTolerance, BinningValue bValue,
                        bool phiCenter = false,
                        std::size_t forceEquidistant = 2u) {
  // Return object (empty for the moment)
  AxisDefs aDef;
  Acts::detail::AxisBoundaryType aType =
      (bValue == binPhi) ? Acts::detail::AxisBoundaryType::Closed
                         : Acts::detail::AxisBoundaryType::Bound;

  std::sort(
      clusters.begin(), clusters.end(),
      [&](const auto& a, const auto& b) { return (a.center < b.center); });

  // Clusterize into diffs
  std::vector<ActsScalar> deltaValues;
  for (auto [ic, c] : enumerate(clusters)) {
    if (ic > 0u) {
      deltaValues.push_back(c.center - clusters[ic - 1u].center);
    }
  }

  // Clusterize the deltas
  auto deltaClusters = clusterize(deltaValues, diffTolerance);

  if (deltaClusters.size() <= forceEquidistant) {
    std::sort(
        deltaClusters.begin(), deltaClusters.end(),
        [&](const auto& a, const auto& b) { return (a.center < b.center); });
    deltaClusters = {deltaClusters[0u]};
  }

  // Nice equidistant binning
  if (deltaClusters.size() == 1u) {
    auto lowest = clusters[0u];
    std::sort(lowest.values.begin(), lowest.values.end());
    auto highest = clusters.back();
    std::sort(highest.values.begin(), highest.values.end());
    ActsScalar min = lowest.values[0u];
    ActsScalar max = highest.values.back();
    std::size_t bins = clusters.size() - 1u;

    if (bValue == binPhi) {
      if (phiCenter) {
        ActsScalar stepPhi = (max - min) / bins;
        min -= 0.5 * stepPhi;
        max += 0.5 * stepPhi;
        bins += 1u;
      }
      // Check for phi for proper closure
      if (std::abs((max - min) - 2 * M_PI) < 0.001) {
        max = min + 2 * M_PI;
      }
    }
    aDef.eqAxisDefs = std::tie(min, max, bins, aType);
  }
  std::vector<ActsScalar> boundaries = {};
  std::for_each(clusters.begin(), clusters.end(),
                [&](const auto& c) { boundaries.push_back(c.center); });
  if (bValue == binPhi) {
    boundaries.push_back(-M_PI);
    boundaries.push_back(M_PI);
  }
  std::sort(boundaries.begin(), boundaries.end());
  aDef.varAxisDefs = std::tie(boundaries, aType);

  return aDef;
}

/// @brief Helper method to fill the indices
///
/// @tparam grid_attacher type, contains grid and casts
///
/// @param gctx the geometry contect
/// @param surfacesPolyhedrons the prepared surfaces with polyhedrons
/// @param [in, out] gAttacher the grid attacher that is filled
///
template <typename grid_attacher>
void fillSurfaces(const GeometryContext& gctx,
                  const std::vector<SurfacPolyhedron>& surfaces,
                  grid_attacher& gAttacher) {
  // Retrive the grid from the attacher
  auto& grid = gAttacher.grid;
  for (auto [is, s] : enumerate(surfaces)) {
    // The vertices for the filling
    std::vector<Vector3> fillVertices = std::get<Polyhedron>(s).vertices;
    const Surface* surface = std::get<const Surface*>(s);
    // Add the surface center
    fillVertices.push_back(surface->center(gctx));
    // Now loop over the fill vertices
    for (const auto& vtx : fillVertices) {
      auto cpos = gAttacher.castPosition(vtx);
      //  Get the reference to the grid entry
      auto& indices = grid.atPosition(cpos);
      // Fill
      if (std::find(indices.begin(), indices.end(), is) == indices.end()) {
        indices.push_back(is);
      }
    }
  }
  // Filled, now complete with neighbor binnings
  populateNeighborhood(grid);
}

/// @brief Fill the navigation state delegate
///
/// @tparam grid_type the type of the grid (is moved into the delegate)
/// @param gctx the geometry context
/// @param nStateUpdator [in,out] the navigation state update
/// @param grid [in] the grid to be filled
/// @param bValues the binning values
/// @param surfacesPolyhedrons the prepared surfaces with polyhedrons
///
template <typename grid_type>
void fillNavigationStateUpdator(
    const GeometryContext& gctx, ManagedNavigationStateUpdator& managedUpdator,
    grid_type&& grid, const std::array<BinningValue, grid_type::DIM>& bValues,
    const std::vector<SurfacPolyhedron>& surfacesPolyhedrons) {
  // Create grid surface and all portals attacher
  auto gSurfaceAttacher =
      Acts::Experimental::detail::GridSurfaceAttacher<grid_type>(
          {std::move(grid), bValues});

  // Now fill the surfaces
  fillSurfaces(gctx, surfacesPolyhedrons, gSurfaceAttacher);

  // Combine to a navigation state updator
  using NavigationUpdatorImpl =
      Acts::Experimental::detail::NavigationStateUpdatorImpl<
          decltype(gSurfaceAttacher), decltype(portalAttacher)>;

  auto navUpdator = std::make_shared<NavigationUpdatorImpl>(
      std::make_tuple(std::move(gSurfaceAttacher), portalAttacher));

  NavigationStateUpdator nStateUpdator;
  nStateUpdator.connect<&NavigationUpdatorImpl::update>(navUpdator.get());
  managedUpdator.delegate = std::move(nStateUpdator);
  managedUpdator.implementation = navUpdator;
}

/// @brief Generate the navigation state upator
///
/// @param gctx the geometry context
/// @param surfaces the surfaces
/// @param bValues are the binValues
///
/// @return
void generateNavigationStateUpdator(
    const GeometryContext& gctx, ManagedNavigationStateUpdator& managedUpdator,
    const std::vector<std::shared_ptr<Surface>>& surfaces,
    const std::array<BinningValue, 2u>& bValues,
    const std::array<ActsScalar, 2u>& clusterTolerance,
    const std::array<ActsScalar, 2u>& equidistantTolerance,
    bool phiCenter = false) noexcept(false) {

  // Make a provide all updator in case there's only one surface
  if (surfaces.size() == 1u){
    managedUpdator = detail::allPortalsAndSurfaces();
    return;
  }

  // Clusterizing
  std::vector<SurfacPolyhedron> surfacesPolyhedrons;
  std::array<std::vector<ActsScalar>, 2u> values;
  surfacesPolyhedrons.reserve(surfaces.size());
  for (auto s : surfaces) {
    // Get the polyhedron
    auto polyhedron = s->polyhedronRepresentation(gctx, 1);
    auto extent = polyhedron.extent();
    surfacesPolyhedrons.push_back(std::make_tuple(s.get(), polyhedron));
    // Get the parameters
    for (auto [ibv, bv] : enumerate(bValues)) {
      if (bv == binPhi and phiCenter) {
        values[ibv].push_back(VectorHelpers::phi(s->center(gctx)));
      } else {
        values[ibv].push_back(extent.min(bv));
        values[ibv].push_back(extent.max(bv));
      }
    }
  }

  std::array<AxisDefs, 2u> axesDefinitions;
  // Now clusterize
  for (auto [ibv, bv] : enumerate(bValues)) {
    auto clusters = clusterize(values[ibv], clusterTolerance[ibv]);
    axesDefinitions[ibv] =
        axisDefinition(clusters, equidistantTolerance[ibv], bv, phiCenter);
  }

  // Create the grid variants
  if (axesDefinitions[0u].eqAxisDefs != std::nullopt and
      axesDefinitions[1u].eqAxisDefs != std::nullopt) {
    // Both are equidistant
    auto [min0, max0, bins0, type0] = axesDefinitions[0u].eqAxisDefs.value();
    auto [min1, max1, bins1, type1] = axesDefinitions[1u].eqAxisDefs.value();
    if (type0 == Acts::detail::AxisBoundaryType::Bound and
        type1 == Acts::detail::AxisBoundaryType::Closed) {
      AxisEqB axis0(min0, max0, bins0);
      AxisEqC axis1(min1, max1, bins1);
      GridEqBEqC grid(std::make_tuple(std::move(axis0), std::move(axis1)));
      fillNavigationStateUpdator(gctx, managedUpdator, std::move(grid), bValues,
                                 surfacesPolyhedrons);
      return;
    } else if (type1 == Acts::detail::AxisBoundaryType::Bound) {
      AxisEqB axis0(min0, max0, bins0);
      AxisEqB axis1(min1, max1, bins1);
      GridEqBEqB grid(std::make_tuple(std::move(axis0), std::move(axis1)));
      fillNavigationStateUpdator(gctx, managedUpdator, std::move(grid), bValues,
                                 surfacesPolyhedrons);
      return;
    }
  } else if (axesDefinitions[0u].eqAxisDefs != std::nullopt) {
    auto [min0, max0, bins0, type0] = axesDefinitions[0u].eqAxisDefs.value();
    auto [values1, vtype1] = axesDefinitions[1u].varAxisDefs.value();
    if (type0 == Acts::detail::AxisBoundaryType::Bound and
        vtype1 == Acts::detail::AxisBoundaryType::Closed) {
      AxisEqB axis0(min0, max0, bins0);
      AxisVarC axis1(values1);
      GridEqBVarC grid(std::make_tuple(std::move(axis0), std::move(axis1)));
      fillNavigationStateUpdator(gctx, managedUpdator, std::move(grid), bValues,
                                 surfacesPolyhedrons);
      return;
    } else if (vtype1 == Acts::detail::AxisBoundaryType::Bound) {
      AxisEqB axis0(min0, max0, bins0);
      AxisVarB axis1(values1);
      GridEqBVarB grid(std::make_tuple(std::move(axis0), std::move(axis1)));
      fillNavigationStateUpdator(gctx, managedUpdator, std::move(grid), bValues,
                                 surfacesPolyhedrons);
      return;
    }
  } else {
    auto [values0, vtype0] = axesDefinitions[0u].varAxisDefs.value();
    auto [values1, vtype1] = axesDefinitions[1u].varAxisDefs.value();
    if (vtype0 == Acts::detail::AxisBoundaryType::Bound and
        vtype1 == Acts::detail::AxisBoundaryType::Closed) {
      AxisVarB axis0(values0);
      AxisVarC axis1(values1);
      GridVarBVarC grid(std::make_tuple(std::move(axis0), std::move(axis1)));
      fillNavigationStateUpdator(gctx, managedUpdator, std::move(grid), bValues,
                                 surfacesPolyhedrons);

      return;
    } else if (vtype1 == Acts::detail::AxisBoundaryType::Bound) {
      AxisVarB axis0(values0);
      AxisVarB axis1(values1);
      GridVarBVarB grid(std::make_tuple(std::move(axis0), std::move(axis1)));
      fillNavigationStateUpdator(gctx, managedUpdator, std::move(grid), bValues,
                                 surfacesPolyhedrons);

      return;
    }
  }
  throw std::invalid_argument(
      "SurfaceGridGenerator: wrong axes type combination, either BB or BC");
}

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts