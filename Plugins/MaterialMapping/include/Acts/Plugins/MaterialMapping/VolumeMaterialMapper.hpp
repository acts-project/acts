// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// VolumeMaterialMapper.hpp, Acts project MaterialPlugins
///////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "Acts/Material/Material.hpp"
#include "Acts/Plugins/MaterialMapping/AccumulatedVolumeMaterial.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

namespace Acts {

/// @brief These functions serve to concatenate material evaluated at arbitrary
/// space points on a provided set of grid points. The procedure creates mean
/// material values at these points and can be used to create an interpolated
/// grid of material.
/// @note Since these functions just associate an arbitrary cloud of space
/// points with the material evaluated at these points to a user defined grid,
/// this procedure is a pure pre-processing of an applicable material getter.
/// This last step requires additional interpolation between the grid points.
namespace VolumeMaterialMapper {
  using RecordedMaterial = std::vector<std::pair<Material, Vector3D>>;
  using EAxis            = detail::EquidistantAxis;
  using Grid2D = detail::Grid<AccumulatedVolumeMaterial, EAxis, EAxis>;
  using Grid3D = detail::Grid<AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
  using MaterialGrid2D = detail::Grid<ActsVectorF<5>, EAxis, EAxis>;
  using MaterialGrid3D = detail::Grid<ActsVectorF<5>, EAxis, EAxis, EAxis>;

  /// @brief Helper method that creates the cache grid for the mapping. This
  /// grid allows the collection of material at a the anchro points.
  ///
  /// @param [in] gridAxis1 Vector of grid points in the arbitrary first
  /// dimension
  /// @param [in] gridAxis2 Vector of grid points in the arbitrary second
  /// dimension
  ///
  /// @return The grid
  Grid2D
  createGrid(std::vector<double> gridAxis1, std::vector<double> gridAxis2)
  {
    // sort the values
    std::sort(gridAxis1.begin(), gridAxis1.end());
    std::sort(gridAxis2.begin(), gridAxis2.end());
    // Get unique values
    gridAxis1.erase(std::unique(gridAxis1.begin(), gridAxis1.end()),
                    gridAxis1.end());
    gridAxis2.erase(std::unique(gridAxis2.begin(), gridAxis2.end()),
                    gridAxis2.end());
    gridAxis1.shrink_to_fit();
    gridAxis2.shrink_to_fit();
    // get the number of bins
    size_t nBinsAxis1 = gridAxis1.size();
    size_t nBinsAxis2 = gridAxis2.size();

    // get the minimum and maximum
    auto minMaxAxis1 = std::minmax_element(gridAxis1.begin(), gridAxis1.end());
    auto minMaxAxis2 = std::minmax_element(gridAxis2.begin(), gridAxis2.end());
    double minAxis1  = *minMaxAxis1.first;
    double minAxis2  = *minMaxAxis2.first;
    double maxAxis1  = *minMaxAxis1.second;
    double maxAxis2  = *minMaxAxis2.second;
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
  /// grid allows the collection of material at a the anchro points.
  ///
  /// @param [in] gridAxis1 Vector of grid points in the arbitrary first
  /// dimension
  /// @param [in] gridAxis2 Vector of grid points in the arbitrary second
  /// dimension
  /// @param [in] gridAxis3 Vector of grid points in the arbitrary third
  /// dimension
  ///
  /// @return The grid
  Grid3D
  createGrid(std::vector<double> gridAxis1,
             std::vector<double> gridAxis2,
             std::vector<double> gridAxis3)
  {
    // sort the values
    std::sort(gridAxis1.begin(), gridAxis1.end());
    std::sort(gridAxis2.begin(), gridAxis2.end());
    std::sort(gridAxis3.begin(), gridAxis3.end());
    // Get unique values
    gridAxis1.erase(std::unique(gridAxis1.begin(), gridAxis1.end()),
                    gridAxis1.end());
    gridAxis2.erase(std::unique(gridAxis2.begin(), gridAxis2.end()),
                    gridAxis2.end());
    gridAxis3.erase(std::unique(gridAxis3.begin(), gridAxis3.end()),
                    gridAxis3.end());
    gridAxis1.shrink_to_fit();
    gridAxis2.shrink_to_fit();
    gridAxis3.shrink_to_fit();
    // get the number of bins
    size_t nBinsAxis1 = gridAxis1.size();
    size_t nBinsAxis2 = gridAxis2.size();
    size_t nBinsAxis3 = gridAxis3.size();

    // get the minimum and maximum
    auto minMaxAxis1 = std::minmax_element(gridAxis1.begin(), gridAxis1.end());
    auto minMaxAxis2 = std::minmax_element(gridAxis2.begin(), gridAxis2.end());
    auto minMaxAxis3 = std::minmax_element(gridAxis3.begin(), gridAxis3.end());
    double minAxis1  = *minMaxAxis1.first;
    double minAxis2  = *minMaxAxis2.first;
    double minAxis3  = *minMaxAxis3.first;
    double maxAxis1  = *minMaxAxis1.second;
    double maxAxis2  = *minMaxAxis2.second;
    double maxAxis3  = *minMaxAxis3.second;
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
  /// @param [in] concatenateToGridPoints Function that assigns the space point
  /// of @p mPoints to the grid points by its local index
  ///
  /// @return The average material grid decomposed into classification numbers
  MaterialGrid2D
  mapMaterialPoints(
      Grid2D&                 grid,
      const RecordedMaterial& mPoints,
      const std::function<Grid2D::index_t(const Vector3D&, const Grid2D&)>&
          concatenateToGridPoints)
  {
    // Walk over each point
    for (const auto& rm : mPoints) {
      // Search for fitting grid point and accumulate
      Grid2D::index_t index = concatenateToGridPoints(rm.second, grid);
      grid.at(index).accumulate(rm.first);
    }

    // Build material grid
    // Re-build the axes
    Grid2D::point_t min   = grid.getMin();
    Grid2D::point_t max   = grid.getMax();
    Grid2D::index_t nBins = grid.getNBins();

    EAxis axis1(min[0], max[0], nBins[0]);
    EAxis axis2(min[1], max[1], nBins[1]);

    // Build the grid and fill it with data
    MaterialGrid2D mGrid(std::make_tuple(axis1, axis2));
    for (size_t index = 0; index < grid.size(); index++) {
      mGrid.at(index)
          = grid.at(index).average().decomposeIntoClassificationNumbers();
    }

    return mGrid;
  }

  /// @brief Concatenate a set of material at arbitrary space points on a set of
  /// grid points and produces a grid containing the averaged material values.
  ///
  /// @param [in] grid The material collecting grid
  /// @param [in] mPoints The set of material at the space points
  /// @param [in] concatenateToGridPoints Function that assigns the space point
  /// of @p mPoints to the grid points by its local index
  ///
  /// @return The average material grid decomposed into classification numbers
  MaterialGrid3D
  mapMaterialPoints(
      Grid3D&                 grid,
      const RecordedMaterial& mPoints,
      const std::function<Grid3D::index_t(const Vector3D&, const Grid3D&)>&
          concatenateToGridPoints)
  {
    // Walk over each point
    for (const auto& rm : mPoints) {
      // Search for fitting grid point and accumulate
      Grid3D::index_t index = concatenateToGridPoints(rm.second, grid);
      grid.at(index).accumulate(rm.first);
    }

    // Build material grid
    // Re-build the axes
    Grid3D::point_t min   = grid.getMin();
    Grid3D::point_t max   = grid.getMax();
    Grid3D::index_t nBins = grid.getNBins();

    EAxis axis1(min[0], max[0], nBins[0]);
    EAxis axis2(min[1], max[1], nBins[1]);
    EAxis axis3(min[2], max[2], nBins[2]);

    // Build the grid and fill it with data
    MaterialGrid3D mGrid(std::make_tuple(axis1, axis2, axis3));
    for (size_t index = 0; index < grid.size(); index++) {
      mGrid.at(index)
          = grid.at(index).average().decomposeIntoClassificationNumbers();
    }
    return mGrid;
  }

  /// @brief This function creates a discrete material map
  ///
  /// @param [in] gridAxis1 Vector of grid points in the arbitrary first
  /// dimension
  /// @param [in] gridAxis2 Vector of grid points in the arbitrary second
  /// dimension
  /// @param [in] mPoints The set of material at the space points
  /// @param [in] concatenateToGridPoints Function that assigns the space point
  /// of @p mPoints to the grid points by its local index
  ///
  /// @return The map
  MaterialGrid2D
  createMaterialGrid(
      std::vector<double>     gridAxis1,
      std::vector<double>     gridAxis2,
      const RecordedMaterial& mPoints,
      const std::function<Grid2D::index_t(const Vector3D&, const Grid2D&)>&
          concatenateToGridPoints)
  {
    Grid2D grid = createGrid(std::move(gridAxis1), std::move(gridAxis2));
    return mapMaterialPoints(grid, mPoints, concatenateToGridPoints);
  }

  /// @brief This function creates a discrete material map
  ///
  /// @param [in] gridAxis1 Vector of grid points in the arbitrary first
  /// dimension
  /// @param [in] gridAxis2 Vector of grid points in the arbitrary second
  /// dimension
  /// @param [in] gridAxis3 Vector of grid points in the arbitrary third
  /// dimension
  /// @param [in] mPoints The set of material at the space points
  /// @param [in] concatenateToGridPoints Function that assigns the space point
  /// of @p mPoints to the grid points by its local index
  ///
  /// @return The map
  MaterialGrid3D
  createMaterialGrid(
      std::vector<double>     gridAxis1,
      std::vector<double>     gridAxis2,
      std::vector<double>     gridAxis3,
      const RecordedMaterial& mPoints,
      const std::function<Grid3D::index_t(const Vector3D&, const Grid3D&)>&
          concatenateToGridPoints)
  {
    Grid3D grid = createGrid(
        std::move(gridAxis1), std::move(gridAxis2), std::move(gridAxis3));
    return mapMaterialPoints(grid, mPoints, concatenateToGridPoints);
  }
}
}
