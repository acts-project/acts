// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// VolumeMaterialMapper.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Plugins/MaterialMapping/VolumeMaterialMapper.hpp"

Acts::VolumeMaterialMapper::Grid2D
Acts::VolumeMaterialMapper::createGrid(std::vector<double> gridAxis1,
                                       std::vector<double> gridAxis2) const
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
  auto   minMaxAxis1 = std::minmax_element(gridAxis1.begin(), gridAxis1.end());
  auto   minMaxAxis2 = std::minmax_element(gridAxis2.begin(), gridAxis2.end());
  double minAxis1    = *minMaxAxis1.first;
  double minAxis2    = *minMaxAxis2.first;
  double maxAxis1    = *minMaxAxis1.second;
  double maxAxis2    = *minMaxAxis2.second;
  // calculate maxima (add one last bin, because bin value always corresponds to
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

Acts::VolumeMaterialMapper::Grid3D
Acts::VolumeMaterialMapper::createGrid(std::vector<double> gridAxis1,
                                       std::vector<double> gridAxis2,
                                       std::vector<double> gridAxis3) const
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
  auto   minMaxAxis1 = std::minmax_element(gridAxis1.begin(), gridAxis1.end());
  auto   minMaxAxis2 = std::minmax_element(gridAxis2.begin(), gridAxis2.end());
  auto   minMaxAxis3 = std::minmax_element(gridAxis3.begin(), gridAxis3.end());
  double minAxis1    = *minMaxAxis1.first;
  double minAxis2    = *minMaxAxis2.first;
  double minAxis3    = *minMaxAxis3.first;
  double maxAxis1    = *minMaxAxis1.second;
  double maxAxis2    = *minMaxAxis2.second;
  double maxAxis3    = *minMaxAxis3.second;
  // calculate maxima (add one last bin, because bin value always corresponds to
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

Acts::VolumeMaterialMapper::MaterialGrid2D
Acts::VolumeMaterialMapper::mapMaterialPoints(
    Grid2D&                 grid,
    const RecordedMaterial& mPoints,
    const std::function<Grid2D::index_t(const Vector3D&, const Grid2D&)>&
        concatenateToGridPoints) const
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
  for (size_t index = 0; index < grid.size(); index++)
    mGrid.at(index)
        = grid.at(index).average().decomposeIntoClassificationNumbers();

  return mGrid;
}

Acts::VolumeMaterialMapper::MaterialGrid3D
Acts::VolumeMaterialMapper::mapMaterialPoints(
    Grid3D&                 grid,
    const RecordedMaterial& mPoints,
    const std::function<Grid3D::index_t(const Vector3D&, const Grid3D&)>&
        concatenateToGridPoints) const
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
  for (size_t index = 0; index < grid.size(); index++)
    mGrid.at(index)
        = grid.at(index).average().decomposeIntoClassificationNumbers();

  return mGrid;
}