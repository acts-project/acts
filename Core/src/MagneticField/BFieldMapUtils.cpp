// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/BFieldMapUtils.hpp"

#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/MagneticField/ToroidalField.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <stdexcept>

namespace Acts {

using VectorHelpers::perp;
using VectorHelpers::phi;

InterpolatedBFieldMap<
    Grid<Vector2, Axis<AxisType::Equidistant>, Axis<AxisType::Equidistant>>>
fieldMapRZ(const std::function<std::size_t(std::array<std::size_t, 2> binsRZ,
                                           std::array<std::size_t, 2> nBinsRZ)>&
               localToGlobalBin,
           std::vector<double> rPos, std::vector<double> zPos,
           const std::vector<Vector2>& bField, double lengthUnit,
           double BFieldUnit, bool firstQuadrant) {
  // [1] Create Grid
  const auto [rMin, rMax, rBinCount] = detail::getMinMaxAndBinCount(rPos);
  auto [zMin, zMax, zBinCount] = detail::getMinMaxAndBinCount(zPos);

  const std::size_t nBinsR = rBinCount;
  std::size_t nBinsZ = zBinCount;

  if (firstQuadrant) {
    zMin = -zPos[nBinsZ - 1];
    nBinsZ = 2 * nBinsZ - 1;
  }

  // Create the axis for the grid
  Axis rAxis(rMin * lengthUnit, rMax * lengthUnit, nBinsR);
  Axis zAxis(zMin * lengthUnit, zMax * lengthUnit, nBinsZ);

  // Create the grid
  Grid grid(Type<Vector2>, std::move(rAxis), std::move(zAxis));
  using Grid_t = decltype(grid);

  // [2] Set the bField values
  const std::array<std::size_t, 2> nIndices = {{rBinCount, zBinCount}};
  for (std::size_t i = 1; i <= nBinsR; ++i) {
    for (std::size_t j = 1; j <= nBinsZ; ++j) {
      Grid_t::index_t indices = {{i, j}};
      // std::vectors begin with 0 and we do not want the user needing to take
      // underflow or overflow bins in account this is why we need to subtract
      // by one
      if (firstQuadrant) {
        std::size_t n = std::abs(static_cast<std::ptrdiff_t>(j) -
                                 static_cast<std::ptrdiff_t>(zBinCount));

        grid.atLocalBins(indices) =
            bField.at(localToGlobalBin({{i - 1, n}}, nIndices)) * BFieldUnit;
      } else {
        grid.atLocalBins(indices) =
            bField.at(localToGlobalBin({{i - 1, j - 1}}, nIndices)) *
            BFieldUnit;
      }
    }
  }
  grid.setExteriorBins(Vector2::Zero());

  // [3] Create the transformation for the position map (x,y,z) -> (r,z)
  auto transformPos = [](const Vector3& pos) {
    return Vector2(perp(pos), pos.z());
  };

  // [4] Create the transformation for the bField map (Br,Bz) -> (Bx,By,Bz)
  auto transformBField = [](const Vector2& field, const Vector3& pos) {
    const double rSinTheta2 = pos.x() * pos.x() + pos.y() * pos.y();
    double cosPhi = 1.;
    double sinPhi = 0.;

    if (rSinTheta2 > std::numeric_limits<double>::min()) {
      const double invRsinTheta = 1. / std::sqrt(rSinTheta2);
      cosPhi = pos.x() * invRsinTheta;
      sinPhi = pos.y() * invRsinTheta;
    }

    return Vector3(field.x() * cosPhi, field.x() * sinPhi, field.y());
  };

  // [5] Create the mapper & BField Service create field mapping
  return InterpolatedBFieldMap<Grid_t>(
      {transformPos, transformBField, std::move(grid)});
}

InterpolatedBFieldMap<
    Grid<Vector3, Axis<AxisType::Equidistant>, Axis<AxisType::Equidistant>,
         Axis<AxisType::Equidistant>>>
fieldMapXYZ(
    const std::function<std::size_t(std::array<std::size_t, 3> binsXYZ,
                                    std::array<std::size_t, 3> nBinsXYZ)>&
        localToGlobalBin,
    std::vector<double> xPos, std::vector<double> yPos,
    std::vector<double> zPos, const std::vector<Vector3>& bField,
    double lengthUnit, double BFieldUnit, bool firstOctant) {
  // [1] Create Grid
  auto [xMin, xMax, xBinCount] = detail::getMinMaxAndBinCount(xPos);
  auto [yMin, yMax, yBinCount] = detail::getMinMaxAndBinCount(yPos);
  auto [zMin, zMax, zBinCount] = detail::getMinMaxAndBinCount(zPos);

  std::size_t nBinsX = xBinCount;
  std::size_t nBinsY = yBinCount;
  std::size_t nBinsZ = zBinCount;

  if (firstOctant) {
    xMin = -xPos[nBinsX - 1];
    nBinsX = 2 * nBinsX - 1;
    yMin = -yPos[nBinsY - 1];
    nBinsY = 2 * nBinsY - 1;
    zMin = -zPos[nBinsZ - 1];
    nBinsZ = 2 * nBinsZ - 1;
  }

  Axis xAxis(xMin * lengthUnit, xMax * lengthUnit, nBinsX);
  Axis yAxis(yMin * lengthUnit, yMax * lengthUnit, nBinsY);
  Axis zAxis(zMin * lengthUnit, zMax * lengthUnit, nBinsZ);
  // Create the grid
  Grid grid(Type<Vector3>, std::move(xAxis), std::move(yAxis),
            std::move(zAxis));
  using Grid_t = decltype(grid);

  // [2] Set the bField values
  const std::array<std::size_t, 3> nIndices = {
      {xBinCount, yBinCount, zBinCount}};

  auto calcAbsDiff = [](std::size_t val, std::size_t binCount) {
    return std::abs(static_cast<std::ptrdiff_t>(val) -
                    static_cast<std::ptrdiff_t>(binCount));
  };

  for (std::size_t i = 1; i <= nBinsX; ++i) {
    for (std::size_t j = 1; j <= nBinsY; ++j) {
      for (std::size_t k = 1; k <= nBinsZ; ++k) {
        Grid_t::index_t indices = {{i, j, k}};
        // std::vectors begin with 0 and we do not want the user needing to take
        // underflow or overflow bins in account this is why we need to subtract
        // by one
        if (firstOctant) {
          const std::size_t l = calcAbsDiff(i, xBinCount);
          const std::size_t m = calcAbsDiff(j, yBinCount);
          const std::size_t n = calcAbsDiff(k, zBinCount);

          grid.atLocalBins(indices) =
              bField.at(localToGlobalBin({{l, m, n}}, nIndices)) * BFieldUnit;
        } else {
          grid.atLocalBins(indices) =
              bField.at(localToGlobalBin({{i - 1, j - 1, k - 1}}, nIndices)) *
              BFieldUnit;
        }
      }
    }
  }
  grid.setExteriorBins(Vector3::Zero());

  // [3] Create the transformation for the position map (x,y,z) -> (r,z)
  auto transformPos = [](const Vector3& pos) { return pos; };

  // [4] Create the transformation for the BField map (Bx,By,Bz) -> (Bx,By,Bz)
  auto transformBField = [](const Vector3& field, const Vector3& /*pos*/) {
    return field;
  };

  // [5] Create the mapper & BField Service create field mapping
  return InterpolatedBFieldMap<Grid_t>(
      {transformPos, transformBField, std::move(grid)});
}

InterpolatedBFieldMap<
    Grid<Vector2, Axis<AxisType::Equidistant>, Axis<AxisType::Equidistant>>>
solenoidFieldMap(const std::pair<double, double>& rLim,
                 const std::pair<double, double>& zLim,
                 const std::pair<std::size_t, std::size_t>& nBins,
                 const SolenoidBField& field) {
  auto [rMin, rMax] = rLim;
  auto [zMin, zMax] = zLim;
  const auto [nBinsR, nBinsZ] = nBins;

  double stepZ = std::abs(zMax - zMin) / (nBinsZ - 1);
  double stepR = std::abs(rMax - rMin) / (nBinsR - 1);
  rMax += stepR;
  zMax += stepZ;

  // Create the axis for the grid
  Axis rAxis(rMin, rMax, nBinsR);
  Axis zAxis(zMin, zMax, nBinsZ);

  // Create the grid
  Grid grid(Type<Vector2>, std::move(rAxis), std::move(zAxis));
  using Grid_t = decltype(grid);

  // Create the transformation for the position map (x,y,z) -> (r,z)
  auto transformPos = [](const Vector3& pos) {
    return Vector2(perp(pos), pos.z());
  };

  // Create the transformation for the bField map (Br,Bz) -> (Bx,By,Bz)
  auto transformBField = [](const Vector2& bField, const Vector3& pos) {
    const double rSinTheta2 = pos.x() * pos.x() + pos.y() * pos.y();
    double cosPhi = 1.;
    double sinPhi = 0.;

    if (rSinTheta2 > std::numeric_limits<double>::min()) {
      const double invRsinTheta = 1. / std::sqrt(rSinTheta2);
      cosPhi = pos.x() * invRsinTheta;
      sinPhi = pos.y() * invRsinTheta;
    }

    return Vector3(bField.x() * cosPhi, bField.x() * sinPhi, bField.y());
  };

  // iterate over all bins, set their value to the solenoid value at their lower
  // left position
  for (std::size_t i = 0; i <= nBinsR + 1; i++) {
    for (std::size_t j = 0; j <= nBinsZ + 1; j++) {
      Grid_t::index_t index({i, j});
      if (i == 0 || j == 0 || i == nBinsR + 1 || j == nBinsZ + 1) {
        // under or overflow bin, set zero
        grid.atLocalBins(index) = Grid_t::value_type(0, 0);
      } else {
        // regular bin, get lower left boundary
        Grid_t::point_t lowerLeft = grid.lowerLeftBinEdge(index);
        // do lookup
        Vector2 B = field.getField(Vector2(lowerLeft[0], lowerLeft[1]));
        grid.atLocalBins(index) = B;
      }
    }
  }

  // Create the mapper & BField Service create field mapping
  InterpolatedBFieldMap<Grid_t> map(
      {transformPos, transformBField, std::move(grid)});
  return map;
}

InterpolatedBFieldMap<
    Grid<Vector3, Axis<AxisType::Equidistant>, Axis<AxisType::Equidistant>,
         Axis<AxisType::Equidistant>>>
toroidalFieldMapCyl(
    const std::pair<double, double>& rLim,
    const std::pair<double, double>& phiLim,
    const std::pair<double, double>& zLim,
    const std::tuple<std::size_t, std::size_t, std::size_t>& nBins,
    const ToroidalField& field) {
  // Create magnetic field context and cache
  MagneticFieldContext ctx;
  auto cache = field.makeCache(ctx);

  auto [rMin, rMax] = rLim;
  auto [pMin, pMax] = phiLim;  // radians
  auto [zMin, zMax] = zLim;
  const auto [nBinsR, nBinsP, nBinsZ] = nBins;

  // Guard against degenerate configs
  if (nBinsR < 2 || nBinsP < 2 || nBinsZ < 2) {
    throw std::invalid_argument(
        "toroidalFieldMapCyl: each dimension needs at least 2 bins");
  }

  // Follow solenoid style: extend upper edge by one step to accommodate
  // overflow bin edges
  double stepR = std::abs(rMax - rMin) / static_cast<double>(nBinsR - 1);
  double stepP = std::abs(pMax - pMin) / static_cast<double>(nBinsP - 1);
  double stepZ = std::abs(zMax - zMin) / static_cast<double>(nBinsZ - 1);
  rMax += stepR;
  pMax += stepP;
  zMax += stepZ;

  // Axes
  Axis rAxis(rMin, rMax, nBinsR);
  Axis pAxis(pMin, pMax, nBinsP);
  Axis zAxis(zMin, zMax, nBinsZ);

  // Grid stores (Br, Bphi, Bz)
  Grid grid(Type<Vector3>, std::move(rAxis), std::move(pAxis),
            std::move(zAxis));
  using Grid_t = decltype(grid);

  // (x,y,z) -> (r,phi,z)
  auto transformPos = [](const Vector3& pos) {
    const double r2 = pos.x() * pos.x() + pos.y() * pos.y();
    const double r = std::sqrt(r2);
    double phi = 0.0;
    if (r2 > std::numeric_limits<double>::min()) {
      phi = std::atan2(pos.y(), pos.x());  // in (-pi, pi]
    }
    return Vector3(r, phi,
                   pos.z());  // note: we reuse Vector3 to carry (r,phi,z)
  };

  // (Br,Bphi,Bz) + Cartesian pos -> (Bx,By,Bz)
  auto transformBField = [](const Vector3& bCyl, const Vector3& pos) {
    const double r2 = pos.x() * pos.x() + pos.y() * pos.y();
    double cosPhi = 1.0, sinPhi = 0.0;
    if (r2 > std::numeric_limits<double>::min()) {
      const double invR = 1.0 / std::sqrt(r2);
      cosPhi = pos.x() * invR;
      sinPhi = pos.y() * invR;
    }
    const double Br = bCyl.x();
    const double Bphi = bCyl.y();
    const double Bz = bCyl.z();
    // Cyl->Cart
    const double Bx = Br * cosPhi - Bphi * sinPhi;
    const double By = Br * sinPhi + Bphi * cosPhi;
    return Vector3(Bx, By, Bz);
  };

  // Fill bins (including under/overflow halo as zeros)
  for (std::size_t ir = 0; ir <= nBinsR + 1; ++ir) {
    for (std::size_t ip = 0; ip <= nBinsP + 1; ++ip) {
      for (std::size_t iz = 0; iz <= nBinsZ + 1; ++iz) {
        Grid_t::index_t index({ir, ip, iz});
        if (ir == 0 || ip == 0 || iz == 0 || ir == nBinsR + 1 ||
            ip == nBinsP + 1 || iz == nBinsZ + 1) {
          grid.atLocalBins(index) = Grid_t::value_type(0.0, 0.0, 0.0);
        } else {
          const Grid_t::point_t ll = grid.lowerLeftBinEdge(index);  // (r,phi,z)
          const double r = ll[0];
          const double phi = ll[1];
          const double z = ll[2];

          const double x = r * std::cos(phi);
          const double y = r * std::sin(phi);

          // Query original field (Cartesian)
          auto res = field.getField(Vector3(x, y, z), cache);

          Vector3 Bxyz(0.0, 0.0, 0.0);
          if (res.ok()) {
            Bxyz = *res;
          }

          // Convert to cylindrical components to store
          double cosPhi =
              (r > std::numeric_limits<double>::min()) ? std::cos(phi) : 1.0;
          double sinPhi =
              (r > std::numeric_limits<double>::min()) ? std::sin(phi) : 0.0;

          const double Br = Bxyz.x() * cosPhi + Bxyz.y() * sinPhi;
          const double Bphi = -Bxyz.x() * sinPhi + Bxyz.y() * cosPhi;
          const double Bz = Bxyz.z();

          grid.atLocalBins(index) = Grid_t::value_type(Br, Bphi, Bz);
        }
      }
    }
  }

  // Mapper
  InterpolatedBFieldMap<Grid_t> map(
      {transformPos, transformBField, std::move(grid)});
  return map;
}

InterpolatedBFieldMap<
    Grid<Vector3, Axis<AxisType::Equidistant>, Axis<AxisType::Equidistant>,
         Axis<AxisType::Equidistant>>>
toroidalFieldMapXYZ(
    const std::pair<double, double>& xLim,
    const std::pair<double, double>& yLim,
    const std::pair<double, double>& zLim,
    const std::tuple<std::size_t, std::size_t, std::size_t>& nBins,
    const ToroidalField& field) {
  // Create magnetic field context and cache
  MagneticFieldContext ctx;
  auto cache = field.makeCache(ctx);

  auto [xMin, xMax] = xLim;
  auto [yMin, yMax] = yLim;
  auto [zMin, zMax] = zLim;
  const auto [nBinsX, nBinsY, nBinsZ] = nBins;

  if (nBinsX < 2 || nBinsY < 2 || nBinsZ < 2) {
    throw std::invalid_argument(
        "toroidalFieldMapXYZ: each dimension needs at least 2 bins");
  }

  double stepX = std::abs(xMax - xMin) / static_cast<double>(nBinsX - 1);
  double stepY = std::abs(yMax - yMin) / static_cast<double>(nBinsY - 1);
  double stepZ = std::abs(zMax - zMin) / static_cast<double>(nBinsZ - 1);
  xMax += stepX;
  yMax += stepY;
  zMax += stepZ;

  Axis xAxis(xMin, xMax, nBinsX);
  Axis yAxis(yMin, yMax, nBinsY);
  Axis zAxis(zMin, zMax, nBinsZ);

  Grid grid(Type<Vector3>, std::move(xAxis), std::move(yAxis),
            std::move(zAxis));
  using Grid_t = decltype(grid);

  // Identity position transform: (x,y,z) -> (x,y,z)
  auto transformPos = [](const Vector3& pos) { return pos; };

  // Identity B-field transform: table already stores (Bx,By,Bz)
  auto transformBField = [](const Vector3& bField, const Vector3& /*pos*/) {
    return bField;
  };

  for (std::size_t ix = 0; ix <= nBinsX + 1; ++ix) {
    for (std::size_t iy = 0; iy <= nBinsY + 1; ++iy) {
      for (std::size_t iz = 0; iz <= nBinsZ + 1; ++iz) {
        Grid_t::index_t index({ix, iy, iz});
        if (ix == 0 || iy == 0 || iz == 0 || ix == nBinsX + 1 ||
            iy == nBinsY + 1 || iz == nBinsZ + 1) {
          grid.atLocalBins(index) = Grid_t::value_type(0.0, 0.0, 0.0);
        } else {
          const Grid_t::point_t ll = grid.lowerLeftBinEdge(index);
          const double x = ll[0];
          const double y = ll[1];
          const double z = ll[2];

          auto res = field.getField(Vector3(x, y, z), cache);
          Vector3 B(0.0, 0.0, 0.0);
          if (res.ok()) {
            B = *res;  // already (Bx,By,Bz)
          }
          grid.atLocalBins(index) = Grid_t::value_type(B[0], B[1], B[2]);
        }
      }
    }
  }

  InterpolatedBFieldMap<Grid_t> map(
      {transformPos, transformBField, std::move(grid)});
  return map;
}

}  // namespace Acts
