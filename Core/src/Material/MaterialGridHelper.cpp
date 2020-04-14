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

Grid2D::index_t mapMaterial2D(const Acts::Vector3D& matPos,
                              const Grid2D& grid) {
  double dist = std::numeric_limits<double>::max();
  size_t indexX = 1, indexY = 1;
  // Loop through the first index
  for (size_t i = 1; i < grid.numLocalBins()[0] + 1; i++) {
    double dX = grid.lowerLeftBinEdge({{i, indexY}})[0] - matPos.x();
    if (abs(dX) < dist) {
      // Store distance and index
      dist = abs(dX);
      indexX = i;
    } else {
      break;
    }
  }
  dist = std::numeric_limits<double>::max();
  // Loop through the second index
  for (size_t i = 1; i < grid.numLocalBins()[1] + 1; i++) {
    double dY = grid.lowerLeftBinEdge({{indexX, i}})[1] - matPos.y();
    if (abs(dY) < dist) {
      // Store distance and index
      dist = abs(dY);
      indexY = i;
    } else {
      break;
    }
  }
  return {{indexX, indexY}};
}

Grid3D::index_t mapMaterial3D(const Acts::Vector3D& matPos,
                              const Grid3D& grid) {
  double dist = std::numeric_limits<double>::max();
  size_t indexX = 1, indexY = 1, indexZ = 1;
  double dX = 0, dY = 0, dZ = 0;
  // Loop through the first index
  for (size_t i = 1; i < grid.numLocalBins()[0] + 1; i++) {
    dX = grid.lowerLeftBinEdge({{i, indexY, indexZ}})[0] - matPos.x();
    if (abs(dX) < dist) {
      // Store distance and index
      dist = abs(dX);
      indexX = i;
    } else {
      break;
    }
  }
  dist = std::numeric_limits<double>::max();
  // Loop through the second index
  for (size_t i = 1; i < grid.numLocalBins()[1] + 1; i++) {
    dY = grid.lowerLeftBinEdge({{indexX, i, indexZ}})[1] - matPos.y();
    if (abs(dY) < dist) {
      // Store distance and index
      dist = abs(dY);
      indexY = i;
    } else {
      break;
    }
  }
  dist = std::numeric_limits<double>::max();
  // Loop through the third index
  for (size_t i = 1; i < grid.numLocalBins()[2] + 1; i++) {
    dZ = grid.lowerLeftBinEdge({{indexX, indexY, i}})[2] - matPos.z();
    if (abs(dZ) < dist) {
      // Store distance and index
      dist = abs(dZ);
      indexZ = i;
    } else {
      break;
    }
  }
  return {{indexX, indexY, indexZ}};
}

Grid3D createGrid(
    const BinUtility& bins,
    std::function<Acts::Vector3D(Acts::Vector3D)>& transfoGlobalToLocal) {
  auto bu = bins.binningData();
  // First we nee to create the 3 axis
  std::array<double, 3> gridAxis1;
  std::array<double, 3> gridAxis2;
  std::array<double, 3> gridAxis3;
  bool iscylinder = false;
  bool iscube = false;
  // loop trought the binned data and create the corresponding axis
  for (size_t b = 0; b < bu.size(); b++) {
    switch (bu[b].binvalue) {
      case binX:
        iscube = true;
        gridAxis1[0] = bu[b].min;
        gridAxis1[1] = bu[b].max;
        gridAxis1[2] = bu[b].bins();
        break;

      case binY:
        iscube = true;
        gridAxis2[0] = bu[b].min;
        gridAxis2[1] = bu[b].max;
        gridAxis2[2] = bu[b].bins();
        break;

      case binR:
        iscylinder = true;
        gridAxis1[0] = bu[b].min;
        gridAxis1[1] = bu[b].max;
        gridAxis1[2] = bu[b].bins();
        break;

      case binPhi:
        iscylinder = true;
        gridAxis2[0] = bu[b].min;
        gridAxis2[1] = bu[b].max;
        gridAxis2[2] = bu[b].bins();
        break;

      case binZ:
        gridAxis3[0] = bu[b].min;
        gridAxis3[1] = bu[b].max;
        gridAxis3[2] = bu[b].bins();
        break;

      case binRPhi:
      case binEta:
      case binH:
      case binMag:
      case binValues:
        throw std::invalid_argument(
            "Incorrect bin, should be x,y,z or r,phi,z");
        break;
    }
  }
  if (!(iscylinder || iscube) || (iscylinder && iscube)) {
    throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
  }
  // create the appropriate Global to local transform
  if (iscylinder) {
    transfoGlobalToLocal = [](Acts::Vector3D pos) -> Acts::Vector3D {
      return {VectorHelpers::perp(pos), VectorHelpers::phi(pos), pos.z()};
    };
  }
  if (iscube) {
    transfoGlobalToLocal = [](Acts::Vector3D pos) -> Acts::Vector3D {
      return {pos.x(), pos.y(), pos.z()};
    };
  }
  // return the grid
  return (Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2),
                           std::move(gridAxis3)));
}

MaterialGrid2D mapMaterialPoints(
    Grid2D& grid, const RecordedMaterialPoint& mPoints,
    std::function<Acts::Vector3D(Acts::Vector3D)>& transfoGlobalToLocal,
    const std::function<Grid2D::index_t(const Acts::Vector3D&, const Grid2D&)>&
        matchToGridPoint) {
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    Grid2D::index_t index =
        matchToGridPoint(transfoGlobalToLocal(rm.second), grid);
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
    std::function<Acts::Vector3D(Acts::Vector3D)>& transfoGlobalToLocal,
    const std::function<Grid3D::index_t(const Acts::Vector3D&, const Grid3D&)>&
        matchToGridPoint) {
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    Grid3D::index_t index =
        matchToGridPoint(transfoGlobalToLocal(rm.second), grid);
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
