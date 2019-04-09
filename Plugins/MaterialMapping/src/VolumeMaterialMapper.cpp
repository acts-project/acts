// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/MaterialMapping/VolumeMaterialMapper.hpp"
#include <iostream>
#include <limits>
#include "Acts/Plugins/MaterialMapping/AccumulatedVolumeMaterial.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/MaterialMapUtils.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace {
using RecordedMaterial = std::vector<std::pair<Acts::Material, Acts::Vector3D>>;
using EAxis            = Acts::detail::EquidistantAxis;
using Grid2D
    = Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis>;
using Grid3D
    = Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
using MaterialGrid2D = Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis>;
using MaterialGrid3D
    = Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis, EAxis>;

/// @brief Helper method that creates the cache grid for the mapping. This
/// grid allows the collection of material at a the anchor points.
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
///
/// @return The grid
Grid2D
createGrid(std::array<double, 3> gridAxis1, std::array<double, 3> gridAxis2)
{
  // get the number of bins
  size_t nBinsAxis1 = gridAxis1[2];
  size_t nBinsAxis2 = gridAxis2[2];

  // get the minimum and maximum
  double minAxis1 = gridAxis1[0];
  double minAxis2 = gridAxis2[0];
  double maxAxis1 = gridAxis1[1];
  double maxAxis2 = gridAxis2[1];
  // calculate maxima (add one last bin, because bin value always corresponds
  // to
  // left boundary)
  double stepAxis1 = std::fabs(maxAxis1 - minAxis1) / (nBinsAxis1 - 1);
  double stepAxis2 = std::fabs(maxAxis2 - minAxis2) / (nBinsAxis2 - 1);
  maxAxis1 += stepAxis1;
  maxAxis2 += stepAxis2;

  // Create the axis for the grid
  EAxis axis1(minAxis1, maxAxis1, nBinsAxis1);
  EAxis axis2(minAxis2, maxAxis2, nBinsAxis2);

  // The material mapping grid
  return Grid2D(std::make_tuple(std::move(axis1), std::move(axis2)));
}

/// @brief Helper method that creates the cache grid for the mapping. This
/// grid allows the collection of material at a the anchor points.
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @param [in] gridAxis3 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
///
/// @return The grid
Grid3D
createGrid(std::array<double, 3> gridAxis1,
           std::array<double, 3> gridAxis2,
           std::array<double, 3> gridAxis3)
{
  // get the number of bins
  size_t nBinsAxis1 = gridAxis1[2];
  size_t nBinsAxis2 = gridAxis2[2];
  size_t nBinsAxis3 = gridAxis3[2];

  // get the minimum and maximum
  double minAxis1 = gridAxis1[0];
  double minAxis2 = gridAxis2[0];
  double minAxis3 = gridAxis3[0];
  double maxAxis1 = gridAxis1[1];
  double maxAxis2 = gridAxis2[1];
  double maxAxis3 = gridAxis3[1];
  // calculate maxima (add one last bin, because bin value always corresponds
  // to
  // left boundary)
  double stepAxis1 = std::fabs(maxAxis1 - minAxis1) / (nBinsAxis1 - 1);
  double stepAxis2 = std::fabs(maxAxis2 - minAxis2) / (nBinsAxis2 - 1);
  double stepAxis3 = std::fabs(maxAxis3 - minAxis3) / (nBinsAxis3 - 1);
  maxAxis1 += stepAxis1;
  maxAxis2 += stepAxis2;
  maxAxis3 += stepAxis3;

  // Create the axis for the grid
  EAxis axis1(minAxis1, maxAxis1, nBinsAxis1);
  EAxis axis2(minAxis2, maxAxis2, nBinsAxis2);
  EAxis axis3(minAxis3, maxAxis3, nBinsAxis3);

  // The material mapping grid
  return Grid3D(
      std::make_tuple(std::move(axis1), std::move(axis2), std::move(axis3)));
}

/// @brief Concatenate a set of material at arbitrary space points on a set of
/// grid points and produces a grid containing the averaged material values.
///
/// @param [in] grid The material collecting grid
/// @param [in] mPoints The set of material at the space points
/// @param [in] matchToGridPoint Function that assigns the space point
/// of @p mPoints to the grid points by its local index
///
/// @return The average material grid decomposed into classification numbers
MaterialGrid2D
mapMaterialPoints(
    Grid2D&                 grid,
    const RecordedMaterial& mPoints,
    const std::function<Grid2D::index_t(const Acts::Vector3D&, const Grid2D&)>&
        matchToGridPoint)
{
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    Grid2D::index_t index = matchToGridPoint(rm.second, grid);
    grid.atLocalBins(index).accumulate(rm.first);
  }

  // Build material grid
  // Re-build the axes
  Grid2D::point_t min   = grid.minPosition();
  Grid2D::point_t max   = grid.maxPosition();
  Grid2D::index_t nBins = grid.numLocalBins();

  EAxis axis1(min[0], max[0], nBins[0]);
  EAxis axis2(min[1], max[1], nBins[1]);

  // Build the grid and fill it with data
  MaterialGrid2D mGrid(std::make_tuple(axis1, axis2));
  for (size_t index = 0; index < grid.size(); index++) {
    mGrid.at(index) = grid.at(index).average().classificationNumbers();
  }

  return mGrid;
}

/// @brief Concatenate a set of material at arbitrary space points on a set of
/// grid points and produces a grid containing the averaged material values.
///
/// @param [in] grid The material collecting grid
/// @param [in] mPoints The set of material at the space points
/// @param [in] matchToGridPoint Function that assigns the space point
/// of @p mPoints to the grid points by its local index
///
/// @return The average material grid decomposed into classification numbers
MaterialGrid3D
mapMaterialPoints(
    Grid3D&                 grid,
    const RecordedMaterial& mPoints,
    const std::function<Grid3D::index_t(const Acts::Vector3D&, const Grid3D&)>&
        matchToGridPoint)
{
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    Grid3D::index_t index = matchToGridPoint(rm.second, grid);
    grid.atLocalBins(index).accumulate(rm.first);
  }

  // Build material grid
  // Re-build the axes
  Grid3D::point_t min   = grid.minPosition();
  Grid3D::point_t max   = grid.maxPosition();
  Grid3D::index_t nBins = grid.numLocalBins();

  EAxis axis1(min[0], max[0], nBins[0]);
  EAxis axis2(min[1], max[1], nBins[1]);
  EAxis axis3(min[2], max[2], nBins[2]);

  // Build the grid and fill it with data
  MaterialGrid3D mGrid(std::make_tuple(axis1, axis2, axis3));
  for (size_t index = 0; index < grid.size(); index++) {
    mGrid.at(index) = grid.at(index).average().classificationNumbers();
  }
  return mGrid;
}
}  // namespace

MaterialGrid2D
Acts::createMaterialGrid(
    std::array<double, 3> gridAxis1,
    std::array<double, 3> gridAxis2,
    const RecordedMaterial& mPoints,
    const std::function<Grid2D::index_t(const Vector3D&, const Grid2D&)>&
        matchToGridPoint)
{
  Grid2D grid = createGrid(std::move(gridAxis1), std::move(gridAxis2));
  return mapMaterialPoints(grid, mPoints, matchToGridPoint);
}

MaterialGrid3D
Acts::createMaterialGrid(
    std::array<double, 3> gridAxis1,
    std::array<double, 3> gridAxis2,
    std::array<double, 3> gridAxis3,
    const RecordedMaterial& mPoints,
    const std::function<Grid3D::index_t(const Vector3D&, const Grid3D&)>&
        matchToGridPoint)
{
  Grid3D grid = createGrid(
      std::move(gridAxis1), std::move(gridAxis2), std::move(gridAxis3));
  return mapMaterialPoints(grid, mPoints, matchToGridPoint);
}