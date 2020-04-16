// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialGridHelper.hpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Material/MaterialGridHelper.hpp"

namespace Acts {

Grid2D createGrid(std::array<double, 3> gridAxis1,
                  std::array<double, 3> gridAxis2) {
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

Grid3D createGrid(std::array<double, 3> gridAxis1,
                  std::array<double, 3> gridAxis2,
                  std::array<double, 3> gridAxis3) {
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
Grid2D createGrid2D(
    const BinUtility& bins,
    std::function<Acts::Vector2D(Acts::Vector3D)>& transfoGlobalToLocal) {
  auto bu = bins.binningData();
  // First we nee to create the 2 axis
  std::array<double, 3> gridAxis1;
  std::array<double, 3> gridAxis2;

  std::vector<Acts::BinningValue> bintype;

  // loop trought the binned data and create the determine the bin type
  for (size_t b = 0; b < bu.size(); b++) {
    if (bu[b].binvalue == binRPhi || bu[b].binvalue == binEta ||
        bu[b].binvalue == binH || bu[b].binvalue == binMag ||
        bu[b].binvalue == binValues) {
      throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
      break;
    }
    bintype.push_back(bu[b].binvalue);
  }

  // Create a 2D Grid in XY
  if (std::find(bintype.begin(), bintype.end(), binX) != bintype.end() &&
      std::find(bintype.begin(), bintype.end(), binY) != bintype.end()) {
    int index1 =
        std::find(bintype.begin(), bintype.end(), binX) - bintype.begin();
    int index2 =
        std::find(bintype.begin(), bintype.end(), binY) - bintype.begin();

    gridAxis1[0] = bu[index1].min;
    gridAxis1[1] = bu[index1].max;
    gridAxis1[2] = bu[index1].bins();

    gridAxis2[0] = bu[index2].min;
    gridAxis2[1] = bu[index2].max;
    gridAxis2[2] = bu[index2].bins();

    transfoGlobalToLocal = [](Acts::Vector3D pos) -> Acts::Vector2D {
      return {pos.x(), pos.y()};
    };
    return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2)));
  }

  // Create a 2D Grid in XZ
  if (std::find(bintype.begin(), bintype.end(), binX) != bintype.end() &&
      std::find(bintype.begin(), bintype.end(), binZ) != bintype.end()) {
    int index1 =
        std::find(bintype.begin(), bintype.end(), binX) - bintype.begin();
    int index2 =
        std::find(bintype.begin(), bintype.end(), binZ) - bintype.begin();

    gridAxis1[0] = bu[index1].min;
    gridAxis1[1] = bu[index1].max;
    gridAxis1[2] = bu[index1].bins();

    gridAxis2[0] = bu[index2].min;
    gridAxis2[1] = bu[index2].max;
    gridAxis2[2] = bu[index2].bins();

    transfoGlobalToLocal = [](Acts::Vector3D pos) -> Acts::Vector2D {
      return {pos.x(), pos.z()};
    };
    return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2)));
  }

  // Create a 2D Grid in YZ
  if (std::find(bintype.begin(), bintype.end(), binY) != bintype.end() &&
      std::find(bintype.begin(), bintype.end(), binZ) != bintype.end()) {
    int index1 =
        std::find(bintype.begin(), bintype.end(), binY) - bintype.begin();
    int index2 =
        std::find(bintype.begin(), bintype.end(), binZ) - bintype.begin();

    gridAxis1[0] = bu[index1].min;
    gridAxis1[1] = bu[index1].max;
    gridAxis1[2] = bu[index1].bins();

    gridAxis2[0] = bu[index2].min;
    gridAxis2[1] = bu[index2].max;
    gridAxis2[2] = bu[index2].bins();

    transfoGlobalToLocal = [](Acts::Vector3D pos) -> Acts::Vector2D {
      return {pos.y(), pos.z()};
    };
    return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2)));
  }

  // Create a 2D Grid in RPhi
  if (std::find(bintype.begin(), bintype.end(), binR) != bintype.end() &&
      std::find(bintype.begin(), bintype.end(), binPhi) != bintype.end()) {
    int index1 =
        std::find(bintype.begin(), bintype.end(), binR) - bintype.begin();
    int index2 =
        std::find(bintype.begin(), bintype.end(), binPhi) - bintype.begin();

    gridAxis1[0] = bu[index1].min;
    gridAxis1[1] = bu[index1].max;
    gridAxis1[2] = bu[index1].bins();

    gridAxis2[0] = bu[index2].min;
    gridAxis2[1] = bu[index2].max;
    gridAxis2[2] = bu[index2].bins();

    transfoGlobalToLocal = [](Acts::Vector3D pos) -> Acts::Vector2D {
      return {VectorHelpers::perp(pos), VectorHelpers::phi(pos)};
    };
    return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2)));
  }

  // Create a 2D Grid in RZ
  if (std::find(bintype.begin(), bintype.end(), binR) != bintype.end() &&
      std::find(bintype.begin(), bintype.end(), binZ) != bintype.end()) {
    int index1 =
        std::find(bintype.begin(), bintype.end(), binR) - bintype.begin();
    int index2 =
        std::find(bintype.begin(), bintype.end(), binZ) - bintype.begin();

    gridAxis1[0] = bu[index1].min;
    gridAxis1[1] = bu[index1].max;
    gridAxis1[2] = bu[index1].bins();

    gridAxis2[0] = bu[index2].min;
    gridAxis2[1] = bu[index2].max;
    gridAxis2[2] = bu[index2].bins();

    transfoGlobalToLocal = [](Acts::Vector3D pos) -> Acts::Vector2D {
      return {VectorHelpers::perp(pos), pos.z()};
    };
    return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2)));
  }

  // Create a 2D Grid in PhiZ
  if (std::find(bintype.begin(), bintype.end(), binPhi) != bintype.end() &&
      std::find(bintype.begin(), bintype.end(), binZ) != bintype.end()) {
    int index1 =
        std::find(bintype.begin(), bintype.end(), binPhi) - bintype.begin();
    int index2 =
        std::find(bintype.begin(), bintype.end(), binZ) - bintype.begin();

    gridAxis1[0] = bu[index1].min;
    gridAxis1[1] = bu[index1].max;
    gridAxis1[2] = bu[index1].bins();

    gridAxis2[0] = bu[index2].min;
    gridAxis2[1] = bu[index2].max;
    gridAxis2[2] = bu[index2].bins();

    transfoGlobalToLocal = [](Acts::Vector3D pos) -> Acts::Vector2D {
      return {VectorHelpers::phi(pos), pos.z()};
    };
    return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2)));
  }

  throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
  return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2)));
}

Grid3D createGrid3D(
    const BinUtility& bins,
    std::function<Acts::Vector3D(Acts::Vector3D)>& transfoGlobalToLocal) {
  auto bu = bins.binningData();
  // First we nee to create the 3 axis
  std::array<double, 3> gridAxis1;
  std::array<double, 3> gridAxis2;
  std::array<double, 3> gridAxis3;

  std::vector<Acts::BinningValue> bintype;

  // loop trought the binned data and create the determine the bin type
  for (size_t b = 0; b < bu.size(); b++) {
    if (bu[b].binvalue == binRPhi || bu[b].binvalue == binEta ||
        bu[b].binvalue == binH || bu[b].binvalue == binMag ||
        bu[b].binvalue == binValues) {
      throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
      break;
    }
    bintype.push_back(bu[b].binvalue);
  }

  // Create a 3D Grid in XYZ
  if (std::find(bintype.begin(), bintype.end(), binX) != bintype.end() &&
      std::find(bintype.begin(), bintype.end(), binY) != bintype.end() &&
      std::find(bintype.begin(), bintype.end(), binZ) != bintype.end()) {
    int index1 =
        std::find(bintype.begin(), bintype.end(), binX) - bintype.begin();
    int index2 =
        std::find(bintype.begin(), bintype.end(), binY) - bintype.begin();
    int index3 =
        std::find(bintype.begin(), bintype.end(), binZ) - bintype.begin();

    gridAxis1[0] = bu[index1].min;
    gridAxis1[1] = bu[index1].max;
    gridAxis1[2] = bu[index1].bins();

    gridAxis2[0] = bu[index2].min;
    gridAxis2[1] = bu[index2].max;
    gridAxis2[2] = bu[index2].bins();

    gridAxis3[0] = bu[index3].min;
    gridAxis3[1] = bu[index3].max;
    gridAxis3[2] = bu[index3].bins();

    transfoGlobalToLocal = [](Acts::Vector3D pos) -> Acts::Vector3D {
      return {pos.x(), pos.y(), pos.z()};
    };
    return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2),
                             std::move(gridAxis3)));
  }

  // Create a 3D Grid in RPhiZ
  if (std::find(bintype.begin(), bintype.end(), binR) != bintype.end() &&
      std::find(bintype.begin(), bintype.end(), binPhi) != bintype.end()) {
    int index1 =
        std::find(bintype.begin(), bintype.end(), binR) - bintype.begin();
    int index2 =
        std::find(bintype.begin(), bintype.end(), binPhi) - bintype.begin();
    int index3 =
        std::find(bintype.begin(), bintype.end(), binZ) - bintype.begin();

    gridAxis1[0] = bu[index1].min;
    gridAxis1[1] = bu[index1].max;
    gridAxis1[2] = bu[index1].bins();

    gridAxis2[0] = bu[index2].min;
    gridAxis2[1] = bu[index2].max;
    gridAxis2[2] = bu[index2].bins();

    gridAxis3[0] = bu[index3].min;
    gridAxis3[1] = bu[index3].max;
    gridAxis3[2] = bu[index3].bins();

    transfoGlobalToLocal = [](Acts::Vector3D pos) -> Acts::Vector3D {
      return {VectorHelpers::perp(pos), VectorHelpers::phi(pos), pos.z()};
    };
    return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2),
                             std::move(gridAxis3)));
  }

  throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
  return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2),
                           std::move(gridAxis3)));
}

MaterialGrid2D mapMaterialPoints(
    Grid2D& grid, const RecordedMaterialPoint& mPoints,
    std::function<Acts::Vector2D(Acts::Vector3D)>& transfoGlobalToLocal) {
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    Grid2D::index_t index =
        grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(rm.second));
    grid.atLocalBins(index).accumulate(rm.first);
  }

  // Build material grid
  // Re-build the axes
  Grid2D::point_t min = grid.minPosition();
  Grid2D::point_t max = grid.maxPosition();
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

MaterialGrid3D mapMaterialPoints(
    Grid3D& grid, const RecordedMaterialPoint& mPoints,
    std::function<Acts::Vector3D(Acts::Vector3D)>& transfoGlobalToLocal) {
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    Grid3D::index_t index =
        grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(rm.second));
    grid.atLocalBins(index).accumulate(rm.first);
  }

  // Build material grid
  // Re-build the axes
  Grid3D::point_t min = grid.minPosition();
  Grid3D::point_t max = grid.maxPosition();
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

}  // namespace Acts
