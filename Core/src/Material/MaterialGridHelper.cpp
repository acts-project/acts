// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/MaterialGridHelper.hpp"

#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cmath>
#include <stdexcept>
#include <tuple>

Acts::Grid2D Acts::createGrid(std::array<double, 3> gridAxis1,
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
  Acts::EAxis axis1(minAxis1, maxAxis1, nBinsAxis1);
  Acts::EAxis axis2(minAxis2, maxAxis2, nBinsAxis2);

  // The material mapping grid
  return Acts::Grid2D(std::make_tuple(std::move(axis1), std::move(axis2)));
}

Acts::Grid3D Acts::createGrid(std::array<double, 3> gridAxis1,
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
  double stepAxis1 =
      std::fabs(maxAxis1 - minAxis1) / std::max(nBinsAxis1 - 1, size_t(1));
  double stepAxis2 =
      std::fabs(maxAxis2 - minAxis2) / std::max(nBinsAxis2 - 1, size_t(1));
  double stepAxis3 =
      std::fabs(maxAxis3 - minAxis3) / std::max(nBinsAxis3 - 1, size_t(1));
  maxAxis1 += stepAxis1;
  maxAxis2 += stepAxis2;
  maxAxis3 += stepAxis3;

  // Create the axis for the grid
  Acts::EAxis axis1(minAxis1, maxAxis1, nBinsAxis1);
  Acts::EAxis axis2(minAxis2, maxAxis2, nBinsAxis2);
  Acts::EAxis axis3(minAxis3, maxAxis3, nBinsAxis3);

  // The material mapping grid
  return Acts::Grid3D(
      std::make_tuple(std::move(axis1), std::move(axis2), std::move(axis3)));
}

std::function<double(Acts::Vector3)> Acts::globalToLocalFromBin(
    Acts::BinningValue& type) {
  std::function<double(Acts::Vector3)> transfoGlobalToLocal;

  switch (type) {
    case Acts::binX:
      transfoGlobalToLocal = [](const Acts::Vector3& pos) -> double {
        return (pos.x());
      };
      break;

    case Acts::binY:
      transfoGlobalToLocal = [](const Acts::Vector3& pos) -> double {
        return (pos.y());
      };
      break;

    case Acts::binR:
      transfoGlobalToLocal = [](const Acts::Vector3& pos) -> double {
        return (Acts::VectorHelpers::perp(pos));
      };
      break;

    case Acts::binPhi:
      transfoGlobalToLocal = [](const Acts::Vector3& pos) -> double {
        return (Acts::VectorHelpers::phi(pos));
      };
      break;

    case Acts::binZ:
      transfoGlobalToLocal = [](const Acts::Vector3& pos) -> double {
        return (pos.z());
      };
      break;

    case Acts::binRPhi:
    case Acts::binEta:
    case Acts::binH:
    case Acts::binMag:
    case Acts::binValues:
      throw std::invalid_argument("Incorrect bin, should be x,y,z,r,phi,z");
      break;
  }
  return transfoGlobalToLocal;
}

Acts::Grid2D Acts::createGrid2D(
    const Acts::BinUtility& bins,
    std::function<Acts::Vector2(Acts::Vector3)>& transfoGlobalToLocal) {
  auto bu = bins.binningData();
  // First we nee to create the 2 axis
  std::array<double, 3> gridAxis1{};
  std::array<double, 3> gridAxis2{};

  bool isCartesian = false;
  bool isCylindrical = false;

  for (size_t b = 0; b < bu.size(); b++) {
    if (bu[b].binvalue == Acts::binX || bu[b].binvalue == Acts::binY) {
      isCartesian = true;
    }
    if (bu[b].binvalue == Acts::binR || bu[b].binvalue == Acts::binPhi) {
      isCylindrical = true;
    }
  }
  if (!(isCartesian || isCylindrical) || (isCylindrical && isCartesian)) {
    throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
  }

  gridAxis1[0] = bu[0].min;
  gridAxis1[1] = bu[0].max;
  gridAxis1[2] = bu[0].bins();

  gridAxis2[0] = bu[1].min;
  gridAxis2[1] = bu[1].max;
  gridAxis2[2] = bu[1].bins();

  std::function<double(Acts::Vector3)> coord1 =
      globalToLocalFromBin(bu[0].binvalue);
  std::function<double(Acts::Vector3)> coord2 =
      globalToLocalFromBin(bu[1].binvalue);
  Transform3 transfo = bins.transform().inverse();
  transfoGlobalToLocal = [coord1, coord2,
                          transfo](Acts::Vector3 pos) -> Acts::Vector2 {
    pos = transfo * pos;
    return {coord1(pos), coord2(pos)};
  };
  return (Acts::createGrid(gridAxis1, gridAxis2));
}

Acts::Grid3D Acts::createGrid3D(
    const Acts::BinUtility& bins,
    std::function<Acts::Vector3(Acts::Vector3)>& transfoGlobalToLocal) {
  auto bu = bins.binningData();
  // First we nee to create the 3 axis
  std::array<double, 3> gridAxis1{};
  std::array<double, 3> gridAxis2{};
  std::array<double, 3> gridAxis3{};

  bool isCartesian = false;
  bool isCylindrical = false;

  for (size_t b = 0; b < bu.size(); b++) {
    if (bu[b].binvalue == Acts::binX || bu[b].binvalue == Acts::binY) {
      isCartesian = true;
    }
    if (bu[b].binvalue == Acts::binR || bu[b].binvalue == Acts::binPhi) {
      isCylindrical = true;
    }
  }
  if (!(isCartesian || isCylindrical) || (isCylindrical && isCartesian)) {
    throw std::invalid_argument("Incorrect bin, should be x,y,z or r,phi,z");
  }

  gridAxis1[0] = bu[0].min;
  gridAxis1[1] = bu[0].max;
  gridAxis1[2] = bu[0].bins();

  gridAxis2[0] = bu[1].min;
  gridAxis2[1] = bu[1].max;
  gridAxis2[2] = bu[1].bins();

  gridAxis3[0] = bu[2].min;
  gridAxis3[1] = bu[2].max;
  gridAxis3[2] = bu[2].bins();

  std::function<double(Acts::Vector3)> coord1 =
      globalToLocalFromBin(bu[0].binvalue);
  std::function<double(Acts::Vector3)> coord2 =
      globalToLocalFromBin(bu[1].binvalue);
  std::function<double(Acts::Vector3)> coord3 =
      globalToLocalFromBin(bu[2].binvalue);
  Transform3 transfo = bins.transform().inverse();

  transfoGlobalToLocal = [coord1, coord2, coord3,
                          transfo](Acts::Vector3 pos) -> Acts::Vector3 {
    pos = transfo * pos;
    return {coord1(pos), coord2(pos), coord3(pos)};
  };
  return (Acts::createGrid(gridAxis1, gridAxis2, gridAxis3));
}

Acts::MaterialGrid2D Acts::mapMaterialPoints(Acts::Grid2D& grid) {
  // Build material grid
  // Re-build the axes
  Acts::Grid2D::point_t min = grid.minPosition();
  Acts::Grid2D::point_t max = grid.maxPosition();
  Acts::Grid2D::index_t nBins = grid.numLocalBins();

  Acts::EAxis axis1(min[0], max[0], nBins[0]);
  Acts::EAxis axis2(min[1], max[1], nBins[1]);

  // Fill the material Grid by averaging the material in the 2D grid
  Acts::MaterialGrid2D mGrid(std::make_tuple(axis1, axis2));
  for (size_t index = 0; index < grid.size(); index++) {
    mGrid.at(index) = grid.at(index).average().parameters();
  }

  return mGrid;
}

Acts::MaterialGrid3D Acts::mapMaterialPoints(Acts::Grid3D& grid) {
  // Build material grid
  // Re-build the axes
  Acts::Grid3D::point_t min = grid.minPosition();
  Acts::Grid3D::point_t max = grid.maxPosition();
  Acts::Grid3D::index_t nBins = grid.numLocalBins();

  Acts::EAxis axis1(min[0], max[0], nBins[0]);
  Acts::EAxis axis2(min[1], max[1], nBins[1]);
  Acts::EAxis axis3(min[2], max[2], nBins[2]);

  // Fill the material Grid by averaging the material in the 3D grid
  Acts::MaterialGrid3D mGrid(std::make_tuple(axis1, axis2, axis3));
  for (size_t index = 0; index < grid.size(); index++) {
    mGrid.at(index) = grid.at(index).average().parameters();
  }
  return mGrid;
}
