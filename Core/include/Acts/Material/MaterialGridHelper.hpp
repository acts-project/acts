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

#pragma once

#include <stdexcept>

#include "Acts/Material/AccumulatedVolumeMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

namespace {
using RecordedMaterialPoint = std::vector<std::pair<Acts::MaterialProperties, Acts::Vector3D>>;
using EAxis = Acts::detail::EquidistantAxis;
using Grid2D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis>;
using Grid3D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
using MaterialGrid2D = Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis>;
using MaterialGrid3D =
    Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis, EAxis>;
}

namespace Acts {

/// @brief Searcher for closest point in 2D cartesian coordinate
///
/// @param [in] matPos Position of the material
/// @param [in] grid Grid that is used for the look-up
///
/// @return Local grid point with the closest distance to @p matPos
Grid2D::index_t mapMaterialCartesian2D(const Acts::Vector3D& matPos, const Grid3D& grid) {
  double dist = std::numeric_limits<double>::max();
  size_t indexX = 0, indexY = 0;
  // Loop through all elements
  for (size_t i = 1; i < grid.numLocalBins()[0]+1; i++) {
    for (size_t j = 1; j < grid.numLocalBins()[1]+1; j++) {
      // Search the closest distance - elements are ordered
      double dX = grid.lowerLeftBinEdge({{i, j}})[0] - matPos.x();
      double dY = grid.lowerLeftBinEdge({{i, j}})[1] - matPos.y();
      if (std::sqrt(dX * dX + dY * dY) < dist) {
        // Store distance and index
        dist = std::sqrt(dX * dX + dY * dY);
        indexX = i;
        indexY = j;
      } else {  // Break if distance becomes larger
        break;
      }
    }
  }
  return {{indexX, indexY}};
}

/// @brief Searcher for closest point in 3D cartesian coordinate
///
/// @param [in] matPos Position of the material
/// @param [in] grid Grid that is used for the look-up
///
/// @return Local grid point with the closest distance to @p matPos
Grid3D::index_t mapMaterialCartesian3D(const Acts::Vector3D& matPos, const Grid3D& grid) {
  double dist = std::numeric_limits<double>::max();
  size_t indexX = 0, indexY = 0, indexZ = 0;
  // Loop through all elements
  for (size_t i = 1; i < grid.numLocalBins()[0]+1; i++) {
    for (size_t j = 1; j < grid.numLocalBins()[1]+1; j++) {
      for (size_t k = 1; k < grid.numLocalBins()[2]+1; k++) {
        // Search the closest distance - elements are ordered
        double dX = grid.lowerLeftBinEdge({{i, j, k}})[0] - matPos.x();
        double dY = grid.lowerLeftBinEdge({{i, j, k}})[1] - matPos.y();
        double dZ = grid.lowerLeftBinEdge({{i, j, k}})[2] - matPos.z();

        if (std::sqrt(dX * dX + dY * dY + dZ * dZ) < dist) {
          // Store distance and index
          dist = std::sqrt(dX * dX + dY * dY + dZ * dZ);
          indexX = i;
          indexY = j;
          indexZ = k;
        } else {  // Break if distance becomes larger
          break;
        }
      }
    }
  }
  return {{indexX, indexY, indexZ}};
}

/// @brief Searcher for closest point in polar coordinate
///
/// @param [in] matPos Position of the material
/// @param [in] grid Grid that is used for the look-up
///
/// @return Local grid point with the closest distance to @p matPos
Grid2D::index_t mapMaterialCylinderRPhi(const Acts::Vector3D& matPos, const Grid3D& grid) {
  double dist = std::numeric_limits<double>::max();
  size_t indexX = 0, indexY = 0;
  // Loop through all elements
  for (size_t i = 1; i < grid.numLocalBins()[0]+1; i++) {
    for (size_t j = 1; j < grid.numLocalBins()[1]+1; j++) {
      double X = grid.lowerLeftBinEdge({{i, j}})[0]*cos(grid.lowerLeftBinEdge({{i, j}})[1]);
      double Y = grid.lowerLeftBinEdge({{i, j}})[0]*sin(grid.lowerLeftBinEdge({{i, j}})[1]);
      // Search the closest distance - elements are ordered
      double dX = X - matPos.x();
      double dY = Y - matPos.y();
      if (std::sqrt(dX * dX + dY * dY) < dist) {
        // Store distance and index
        dist = std::sqrt(dX * dX + dY * dY );
        indexX = i;
        indexY = j;
      } else {  // Break if distance becomes larger
        // break;
      }
    }
  }
  return {{indexX, indexY}};
}

/// @brief Searcher for closest point in cylindrical coordinate
///
/// @param [in] matPos Position of the material
/// @param [in] grid Grid that is used for the look-up
///
/// @return Local grid point with the closest distance to @p matPos
Grid3D::index_t mapMaterialCylinder3D(const Acts::Vector3D& matPos, const Grid3D& grid) {
  double dist = std::numeric_limits<double>::max();
  size_t indexX = 0, indexY = 0, indexZ = 0;
  // Loop through all elements
  for (size_t i = 1; i < grid.numLocalBins()[0]+1; i++) {
    for (size_t j = 1; j < grid.numLocalBins()[1]+1; j++) {
      for (size_t k = 1; k < grid.numLocalBins()[2]+1; k++) {
        double X = grid.lowerLeftBinEdge({{i, j, k}})[0]*cos(grid.lowerLeftBinEdge({{i, j, k}})[1]);
        double Y = grid.lowerLeftBinEdge({{i, j, k}})[0]*sin(grid.lowerLeftBinEdge({{i, j, k}})[1]);
        // Search the closest distance - elements are ordered
        double dX = X - matPos.x();
        double dY = Y - matPos.y();
        double dZ = grid.lowerLeftBinEdge({{i, j, k}})[2] - matPos.z();
        if (std::sqrt(dX * dX + dY * dY + dZ * dZ) < dist) {
          // Store distance and index
          dist = std::sqrt(dX * dX + dY * dY + dZ * dZ);
          indexX = i;
          indexY = j;
          indexZ = k;
        } else {  // Break if distance becomes larger
          // break;
        }
      }
    }
  }
  return {{indexX, indexY, indexZ}};
}

/// @brief Helper method that creates the cache grid for the mapping. This
/// grid allows the collection of material at a the anchor points.
///
/// @param [in] gridAxis1 Axis data
/// @param [in] gridAxis2 Axis data
/// @note The data of the axes is given in the std::array as {minimum value,
/// maximum value, number of bins}
///
/// @return The grid
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


Grid3D createGrid(const BinUtility &bins,
                  std::function<Acts::Vector3D(Acts::Vector3D)> &transfoGlobalToLocal,
                  std::function<Grid3D::index_t(const Acts::Vector3D&, const Grid3D&)> &mapMaterial) {

  auto bu = bins.binningData();
  std::array<double, 3> gridAxis1;
  std::array<double, 3> gridAxis2;
  std::array<double, 3> gridAxis3;
  bool iscylinder = false;
  bool iscube = false;
  for(size_t b=0; b<bu.size(); b++){
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
      throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
      break;
    }
  }
  if( !(iscylinder || iscube) || (iscylinder && iscube)){
    throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
  }
  if(iscylinder){
    transfoGlobalToLocal = [](Acts::Vector3D pos)->Acts::Vector3D{return{sqrt(pos.x()*pos.x()+pos.y()*pos.y()),
                                                                         atan2(pos.y(),pos.x()), pos.z()}; };
    mapMaterial = Acts::mapMaterialCylinder3D;
  }
  if(iscube){
    transfoGlobalToLocal = [](Acts::Vector3D pos)->Acts::Vector3D{return{pos.x(), pos.y(), pos.z()};};
    mapMaterial = Acts::mapMaterialCartesian3D;
  }
  return(Acts::createGrid(std::move(gridAxis1), std::move(gridAxis2),
                          std::move(gridAxis3)));
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
MaterialGrid2D mapMaterialPoints(
    Grid2D& grid, const Acts::RecordedMaterialPoint& mPoints,
    const std::function<Grid2D::index_t(const Acts::Vector3D&, const Grid2D&)>&
        matchToGridPoint) {
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    Grid2D::index_t index = matchToGridPoint(rm.second, grid);
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

/// @brief Concatenate a set of material at arbitrary space points on a set of
/// grid points and produces a grid containing the averaged material values.
///
/// @param [in] grid The material collecting grid
/// @param [in] mPoints The set of material at the space points
/// @param [in] matchToGridPoint Function that assigns the space point
/// of @p mPoints to the grid points by its local index
///
/// @return The average material grid decomposed into classification numbers
MaterialGrid3D mapMaterialPoints(
    Grid3D& grid, const Acts::RecordedMaterialPoint& mPoints,
    const std::function<Grid3D::index_t(const Acts::Vector3D&, const Grid3D&)>&
        matchToGridPoint) {
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    Grid3D::index_t index = matchToGridPoint(rm.second, grid);
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

}  // Acts
