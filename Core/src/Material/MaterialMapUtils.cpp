// This file is part of the Acts project.
//
// Copyright (C) 2019-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/MaterialMapUtils.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <set>
#include <tuple>
#include <utility>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace {
auto getMinMaxAndBinCount(std::vector<double>& xPos) {
  // sort the values for unique()
  std::ranges::sort(xPos);

  // get the number of bins over unique values
  const std::size_t nBinsX =
      std::distance(xPos.begin(), std::ranges::unique(xPos).end());

  // get the minimum and maximum
  auto [xMin, xMax] = std::ranges::minmax(xPos);

  // calculate maxima (add one last bin, because bin value always corresponds to
  // left boundary)
  const double stepX = (xMax - xMin) / (nBinsX - 1);
  xMax += stepX;

  // Return all values as a tuple
  return std::make_tuple(xMin, xMax, nBinsX);
}
}  // anonymous namespace

auto Acts::materialMapperRZ(
    const std::function<std::size_t(std::array<std::size_t, 2> binsRZ,
                                    std::array<std::size_t, 2> nBinsRZ)>&
        materialVectorToGridMapper,
    std::vector<double> rPos, std::vector<double> zPos,
    const std::vector<Acts::Material>& material, double lengthUnit)
    -> MaterialMapper<
        Grid<Material::ParametersVector, Axis<Acts::AxisType::Equidistant>,
             Axis<Acts::AxisType::Equidistant>>> {
  // [1] Decompose material
  std::vector<Material::ParametersVector> materialVector;
  materialVector.reserve(material.size());

  for (const Material& mat : material) {
    materialVector.push_back(mat.parameters());
  }

  // [2] Create Grid
  const auto [rMin, rMax, nBinsR] = getMinMaxAndBinCount(rPos);
  const auto [zMin, zMax, nBinsZ] = getMinMaxAndBinCount(zPos);

  // Create the axis for the grid
  Axis rAxis(rMin * lengthUnit, rMax * lengthUnit, nBinsR);
  Axis zAxis(zMin * lengthUnit, zMax * lengthUnit, nBinsZ);

  // Create the grid
  Grid grid(Type<Material::ParametersVector>, std::move(rAxis),
            std::move(zAxis));
  using Grid_t = decltype(grid);

  // [3] Set the material values
  const std::array<std::size_t, 2> nIndices = {{nBinsR, nBinsZ}};
  for (std::size_t i = 1; i <= nBinsR; ++i) {
    for (std::size_t j = 1; j <= nBinsZ; ++j) {
      Grid_t::index_t indices = {{i, j}};
      // std::vectors begin with 0 and we do not want the user needing to take
      // underflow or overflow bins in account this is why we need to subtract
      // by one
      grid.atLocalBins(indices) = materialVector.at(
          materialVectorToGridMapper({{i - 1, j - 1}}, nIndices));
    }
  }
  Material::ParametersVector vec;
  vec << std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
      0., 0., 0.;
  grid.setExteriorBins(vec);

  // [4] Create the transformation for the position map (x,y,z) -> (r,z)
  auto transformPos = [](const Vector3& pos) {
    return Vector2(perp(pos), pos.z());
  };

  // [5] Create the mapper & BField Service create material mapping
  return MaterialMapper(transformPos, std::move(grid));
}

auto Acts::materialMapperXYZ(
    const std::function<std::size_t(std::array<std::size_t, 3> binsXYZ,
                                    std::array<std::size_t, 3> nBinsXYZ)>&
        materialVectorToGridMapper,
    std::vector<double> xPos, std::vector<double> yPos,
    std::vector<double> zPos, const std::vector<Material>& material,
    double lengthUnit)
    -> MaterialMapper<Grid<
        Material::ParametersVector, Axis<Acts::AxisType::Equidistant>,
        Axis<Acts::AxisType::Equidistant>, Axis<Acts::AxisType::Equidistant>>> {
  // [1] Decompose material
  std::vector<Material::ParametersVector> materialVector;
  materialVector.reserve(material.size());

  for (const Material& mat : material) {
    materialVector.push_back(mat.parameters());
  }

  // [2] Create Grid
  const auto [xMin, xMax, nBinsX] = getMinMaxAndBinCount(xPos);
  const auto [yMin, yMax, nBinsY] = getMinMaxAndBinCount(yPos);
  const auto [zMin, zMax, nBinsZ] = getMinMaxAndBinCount(zPos);

  // Create the axis for the grid
  Axis xAxis(xMin * lengthUnit, xMax * lengthUnit, nBinsX);
  Axis yAxis(yMin * lengthUnit, yMax * lengthUnit, nBinsY);
  Axis zAxis(zMin * lengthUnit, zMax * lengthUnit, nBinsZ);
  // Create the grid
  Grid grid(Type<Material::ParametersVector>, std::move(xAxis),
            std::move(yAxis), std::move(zAxis));
  using Grid_t = decltype(grid);

  // [3] Set the bField values
  const std::array<std::size_t, 3> nIndices = {{nBinsX, nBinsY, nBinsZ}};
  for (std::size_t i = 1; i <= nBinsX; ++i) {
    for (std::size_t j = 1; j <= nBinsY; ++j) {
      for (std::size_t k = 1; k <= nBinsZ; ++k) {
        Grid_t::index_t indices = {{i, j, k}};
        // std::vectors begin with 0 and we do not want the user needing to take
        // underflow or overflow bins in account this is why we need to subtract
        // by one
        grid.atLocalBins(indices) = materialVector.at(
            materialVectorToGridMapper({{i - 1, j - 1, k - 1}}, nIndices));
      }
    }
  }
  Material::ParametersVector vec;
  vec << std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
      0., 0., 0.;
  grid.setExteriorBins(vec);

  // [4] Create the transformation for the position map (x,y,z) -> (r,z)
  auto transformPos = [](const Vector3& pos) { return pos; };

  // [5] Create the mapper & BField Service create material mapping
  return MaterialMapper(transformPos, std::move(grid));
}
