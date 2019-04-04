// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/MaterialMapUtils.hpp"
#include <iostream>
#include <limits>
#include "Acts/Utilities/Helpers.hpp"
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
    mGrid.at(index) = grid.at(index).average().classificationNumbers();
  }
  return mGrid;
}
}  // namespace

auto
Acts::materialMapperRZ(
    const std::function<size_t(std::array<size_t, 2> binsRZ,
                               std::array<size_t, 2> nBinsRZ)>&
                                materialVectorToGridMapper,
    std::vector<double>         rPos,
    std::vector<double>         zPos,
    std::vector<Acts::Material> material,
    double                      lengthUnit) -> MaterialMapper<MaterialGrid2D>
{
  // [1] Decompose material
  std::vector<ActsVectorF<5>> materialVector;
  materialVector.reserve(material.size());

  for (Material& mat : material) {
    materialVector.push_back(mat.classificationNumbers());
  }

  // [2] Create Grid
  // sort the values
  std::sort(rPos.begin(), rPos.end());
  std::sort(zPos.begin(), zPos.end());
  // Get unique values
  rPos.erase(std::unique(rPos.begin(), rPos.end()), rPos.end());
  zPos.erase(std::unique(zPos.begin(), zPos.end()), zPos.end());
  rPos.shrink_to_fit();
  zPos.shrink_to_fit();
  // get the number of bins
  size_t nBinsR = rPos.size();
  size_t nBinsZ = zPos.size();

  // get the minimum and maximum
  auto   minMaxR = std::minmax_element(rPos.begin(), rPos.end());
  auto   minMaxZ = std::minmax_element(zPos.begin(), zPos.end());
  double rMin    = *minMaxR.first;
  double zMin    = *minMaxZ.first;
  double rMax    = *minMaxR.second;
  double zMax    = *minMaxZ.second;
  // calculate maxima (add one last bin, because bin value always corresponds to
  // left boundary)
  double stepZ = std::fabs(zMax - zMin) / (nBinsZ - 1);
  double stepR = std::fabs(rMax - rMin) / (nBinsR - 1);
  rMax += stepR;
  zMax += stepZ;

  // Create the axis for the grid
  detail::EquidistantAxis rAxis(rMin * lengthUnit, rMax * lengthUnit, nBinsR);
  detail::EquidistantAxis zAxis(zMin * lengthUnit, zMax * lengthUnit, nBinsZ);

  // Create the grid
  using Grid_t = detail::
      Grid<ActsVectorF<5>, detail::EquidistantAxis, detail::EquidistantAxis>;
  Grid_t grid(std::make_tuple(std::move(rAxis), std::move(zAxis)));

  // [3] Set the material values
  for (size_t i = 1; i <= nBinsR; ++i) {
    for (size_t j = 1; j <= nBinsZ; ++j) {
      std::array<size_t, 2> nIndices = {{rPos.size(), zPos.size()}};
      Grid_t::index_t indices = {{i, j}};
      // std::vectors begin with 0 and we do not want the user needing to
      // take underflow or overflow bins in account this is why we need to
      // subtract by one
      grid.at(indices) = materialVector.at(
          materialVectorToGridMapper({{i - 1, j - 1}}, nIndices));
    }
  }
  ActsVectorF<5> vec;
  vec << std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
      0., 0., 0.;
  grid.setExteriorBins(vec);

  // [4] Create the transformation for the position
  // map (x,y,z) -> (r,z)
  auto transformPos
      = [](const Vector3D& pos) { return Vector2D(perp(pos), pos.z()); };

  // [5] Create the mapper & BField Service
  // create material mapping
  return MaterialMapper<MaterialGrid2D>(transformPos, std::move(grid));
}

auto
Acts::materialMapperXYZ(
    const std::function<size_t(std::array<size_t, 3> binsXYZ,
                               std::array<size_t, 3> nBinsXYZ)>&
                          materialVectorToGridMapper,
    std::vector<double>   xPos,
    std::vector<double>   yPos,
    std::vector<double>   zPos,
    std::vector<Material> material,
    double                lengthUnit) -> MaterialMapper<MaterialGrid3D>
{
  // [1] Decompose material
  std::vector<ActsVectorF<5>> materialVector;
  materialVector.reserve(material.size());

  for (Material& mat : material) {
    materialVector.push_back(mat.classificationNumbers());
  }

  // [2] Create Grid
  // Sort the values
  std::sort(xPos.begin(), xPos.end());
  std::sort(yPos.begin(), yPos.end());
  std::sort(zPos.begin(), zPos.end());
  // Get unique values
  xPos.erase(std::unique(xPos.begin(), xPos.end()), xPos.end());
  yPos.erase(std::unique(yPos.begin(), yPos.end()), yPos.end());
  zPos.erase(std::unique(zPos.begin(), zPos.end()), zPos.end());
  xPos.shrink_to_fit();
  yPos.shrink_to_fit();
  zPos.shrink_to_fit();
  // get the number of bins
  size_t nBinsX = xPos.size();
  size_t nBinsY = yPos.size();
  size_t nBinsZ = zPos.size();

  // get the minimum and maximum
  auto minMaxX = std::minmax_element(xPos.begin(), xPos.end());
  auto minMaxY = std::minmax_element(yPos.begin(), yPos.end());
  auto minMaxZ = std::minmax_element(zPos.begin(), zPos.end());
  // Create the axis for the grid
  // get minima
  double xMin = *minMaxX.first;
  double yMin = *minMaxY.first;
  double zMin = *minMaxZ.first;
  // get maxima
  double xMax = *minMaxX.second;
  double yMax = *minMaxY.second;
  double zMax = *minMaxZ.second;
  // calculate maxima (add one last bin, because bin value always corresponds to
  // left boundary)
  double stepZ = std::fabs(zMax - zMin) / (nBinsZ - 1);
  double stepY = std::fabs(yMax - yMin) / (nBinsY - 1);
  double stepX = std::fabs(xMax - xMin) / (nBinsX - 1);
  xMax += stepX;
  yMax += stepY;
  zMax += stepZ;

  detail::EquidistantAxis xAxis(xMin * lengthUnit, xMax * lengthUnit, nBinsX);
  detail::EquidistantAxis yAxis(yMin * lengthUnit, yMax * lengthUnit, nBinsY);
  detail::EquidistantAxis zAxis(zMin * lengthUnit, zMax * lengthUnit, nBinsZ);
  // Create the grid
  using Grid_t = detail::Grid<ActsVectorF<5>,
                              detail::EquidistantAxis,
                              detail::EquidistantAxis,
                              detail::EquidistantAxis>;
  Grid_t grid(
      std::make_tuple(std::move(xAxis), std::move(yAxis), std::move(zAxis)));

  // [3] Set the bField values
  for (size_t i = 1; i <= nBinsX; ++i) {
    for (size_t j = 1; j <= nBinsY; ++j) {
      for (size_t k = 1; k <= nBinsZ; ++k) {
        Grid_t::index_t indices = {{i, j, k}};
        std::array<size_t, 3> nIndices
            = {{xPos.size(), yPos.size(), zPos.size()}};
        // std::vectors begin with 0 and we do not want the user needing to
        // take underflow or overflow bins in account this is why we need to
        // subtract by one
        grid.at(indices) = materialVector.at(
            materialVectorToGridMapper({{i - 1, j - 1, k - 1}}, nIndices));
      }
    }
  }
  ActsVectorF<5> vec;
  vec << std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
      0., 0., 0.;
  grid.setExteriorBins(vec);

  // [4] Create the transformation for the position
  // map (x,y,z) -> (r,z)
  auto transformPos = [](const Vector3D& pos) { return pos; };

  // [5] Create the mapper & BField Service
  // create material mapping
  return MaterialMapper<MaterialGrid3D>(transformPos, std::move(grid));
}

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