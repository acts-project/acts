// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/MaterialMapUtils.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <tuple>
#include <utility>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

auto Acts::materialMapperRZ(
    const std::function<size_t(std::array<size_t, 2> binsRZ,
                               std::array<size_t, 2> nBinsRZ)>&
        materialVectorToGridMapper,
    std::vector<double> rPos, std::vector<double> zPos,
    const std::vector<Acts::Material>& material, double lengthUnit)
    -> MaterialMapper<
        detail::Grid<Material::ParametersVector, detail::EquidistantAxis,
                     detail::EquidistantAxis>> {
  // [1] Decompose material
  std::vector<Material::ParametersVector> materialVector;
  materialVector.reserve(material.size());

  for (const Material& mat : material) {
    materialVector.push_back(mat.parameters());
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
  auto minMaxR = std::minmax_element(rPos.begin(), rPos.end());
  auto minMaxZ = std::minmax_element(zPos.begin(), zPos.end());
  double rMin = *minMaxR.first;
  double zMin = *minMaxZ.first;
  double rMax = *minMaxR.second;
  double zMax = *minMaxZ.second;
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
  using Grid_t = detail::Grid<Material::ParametersVector,
                              detail::EquidistantAxis, detail::EquidistantAxis>;
  Grid_t grid(std::make_tuple(std::move(rAxis), std::move(zAxis)));

  // [3] Set the material values
  for (size_t i = 1; i <= nBinsR; ++i) {
    for (size_t j = 1; j <= nBinsZ; ++j) {
      std::array<size_t, 2> nIndices = {{rPos.size(), zPos.size()}};
      Grid_t::index_t indices = {{i, j}};
      // std::vectors begin with 0 and we do not want the user needing to
      // take underflow or overflow bins in account this is why we need to
      // subtract by one
      grid.atLocalBins(indices) = materialVector.at(
          materialVectorToGridMapper({{i - 1, j - 1}}, nIndices));
    }
  }
  Material::ParametersVector vec;
  vec << std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
      0., 0., 0.;
  grid.setExteriorBins(vec);

  // [4] Create the transformation for the position
  // map (x,y,z) -> (r,z)
  auto transformPos = [](const Vector3& pos) {
    return Vector2(perp(pos), pos.z());
  };

  // [5] Create the mapper & BField Service
  // create material mapping
  return MaterialMapper<
      detail::Grid<Material::ParametersVector, detail::EquidistantAxis,
                   detail::EquidistantAxis>>(transformPos, std::move(grid));
}

auto Acts::materialMapperXYZ(
    const std::function<size_t(std::array<size_t, 3> binsXYZ,
                               std::array<size_t, 3> nBinsXYZ)>&
        materialVectorToGridMapper,
    std::vector<double> xPos, std::vector<double> yPos,
    std::vector<double> zPos, const std::vector<Material>& material,
    double lengthUnit)
    -> MaterialMapper<
        detail::Grid<Material::ParametersVector, detail::EquidistantAxis,
                     detail::EquidistantAxis, detail::EquidistantAxis>> {
  // [1] Decompose material
  std::vector<Material::ParametersVector> materialVector;
  materialVector.reserve(material.size());

  for (const Material& mat : material) {
    materialVector.push_back(mat.parameters());
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
  using Grid_t =
      detail::Grid<Material::ParametersVector, detail::EquidistantAxis,
                   detail::EquidistantAxis, detail::EquidistantAxis>;
  Grid_t grid(
      std::make_tuple(std::move(xAxis), std::move(yAxis), std::move(zAxis)));

  // [3] Set the bField values
  for (size_t i = 1; i <= nBinsX; ++i) {
    for (size_t j = 1; j <= nBinsY; ++j) {
      for (size_t k = 1; k <= nBinsZ; ++k) {
        Grid_t::index_t indices = {{i, j, k}};
        std::array<size_t, 3> nIndices = {
            {xPos.size(), yPos.size(), zPos.size()}};
        // std::vectors begin with 0 and we do not want the user needing to
        // take underflow or overflow bins in account this is why we need to
        // subtract by one
        grid.atLocalBins(indices) = materialVector.at(
            materialVectorToGridMapper({{i - 1, j - 1, k - 1}}, nIndices));
      }
    }
  }
  Material::ParametersVector vec;
  vec << std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
      0., 0., 0.;
  grid.setExteriorBins(vec);

  // [4] Create the transformation for the position
  // map (x,y,z) -> (r,z)
  auto transformPos = [](const Vector3& pos) { return pos; };

  // [5] Create the mapper & BField Service
  // create material mapping
  return MaterialMapper<
      detail::Grid<Material::ParametersVector, detail::EquidistantAxis,
                   detail::EquidistantAxis, detail::EquidistantAxis>>(
      transformPos, std::move(grid));
}
