// This file is part of the Acts project.
//
// Copyright (C) 2020-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>

namespace Acts {

namespace {

/// @brief Helper to retrieve values of an nDim-dimensional normal
/// distribution
/// @note The constant prefactor (2 * pi)^(- nDim / 2) is discarded
///
/// @param args Coordinates where the Gaussian should be evaluated
/// @note args must be in a coordinate system with origin at the mean
/// values of the Gaussian
/// @param cov Covariance matrix
///
/// @return Multivariate Gaussian evaluated at args
template <unsigned int nDim>
double multivariateGaussian(const ActsVector<nDim>& args,
                            const ActsSquareMatrix<nDim>& cov) {
  double exponent = -0.5 * args.transpose().dot(cov.inverse() * args);
  double gaussianDensity = safeExp(exponent) / std::sqrt(cov.determinant());
  return gaussianDensity;
}

}  // namespace

template <int trkGridSize>
double AdaptiveGridTrackDensity<trkGridSize>::getBinCenter(std::int32_t bin,
                                                           double binExtent) {
  return bin * binExtent;
}

template <int trkGridSize>
std::int32_t AdaptiveGridTrackDensity<trkGridSize>::getBin(double value,
                                                           double binExtent) {
  return static_cast<int>(std::floor(value / binExtent - 0.5) + 1);
}

template <int trkGridSize>
Result<double> AdaptiveGridTrackDensity<trkGridSize>::getMaxZPosition(
    MainDensityMap& densityMap) const {
  if (densityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  int zGridPos = -1;
  if (!m_cfg.useHighestSumZPosition) {
    zGridPos = std::distance(
        densityMap.density.begin(),
        std::max_element(densityMap.density.begin(), densityMap.density.end()));
  } else {
    // Get z position with highest density sum
    // of surrounding bins
    zGridPos = highestDensitySumBin(densityMap);
  }

  std::int32_t zBin = densityMap.zBin[zGridPos];
  // Derive corresponding z value
  int sign = (zBin > 0) ? +1 : -1;
  return (zBin + sign * 0.5f) * m_cfg.binSize;
}

template <int trkGridSize>
Result<typename AdaptiveGridTrackDensity<trkGridSize>::ZPositionAndWidth>
AdaptiveGridTrackDensity<trkGridSize>::getMaxZPositionAndWidth(
    MainDensityMap& densityMap) const {
  // Get z maximum value
  auto maxZRes = getMaxZPosition(densityMap);
  if (not maxZRes.ok()) {
    return maxZRes.error();
  }
  double maxZ = *maxZRes;

  // Get seed width estimate
  auto widthRes = estimateSeedWidth(densityMap, maxZ);
  if (not widthRes.ok()) {
    return widthRes.error();
  }
  double width = *widthRes;

  return ZPositionAndWidth{maxZ, width};
}

template <int trkGridSize>
typename AdaptiveGridTrackDensity<trkGridSize>::TrackDensityMap
AdaptiveGridTrackDensity<trkGridSize>::addTrack(
    const BoundTrackParameters& trk, MainDensityMap& densityMap) const {
  ActsVector<2> impactParams = trk.spatialImpactParameters();
  ActsSquareMatrix<2> cov = trk.spatialImpactParameterCovariance().value();

  double d0 = trk.parameters()[0];
  double z0 = trk.parameters()[1];

  // Calculate offset in d direction to central bin at z-axis
  int dOffset = static_cast<int>(std::floor(d0 / m_cfg.binSize - 0.5) + 1);

  // Check if current track does affect grid density
  // in central bins at z-axis
  if (std::abs(dOffset) > (trkGridSize - 1) / 2.) {
    return {};
  }

  // Calculate bin in z
  int zBin = int(z0 / m_cfg.binSize);

  TrackDensityMap trackDensityMap = createTrackGrid(zBin, impactParams, cov);

  std::vector<std::int32_t> zBinValues;

  std::int32_t startEnd = std::int32_t(trkGridSize - 1) / 2;

  for (std::int32_t i = 0; i < trkGridSize; ++i) {
    zBinValues.push_back(static_cast<std::int32_t>(zBin + (i - startEnd)));
  }

  for (std::int32_t i = 0; i < trkGridSize; ++i) {
    std::int32_t z = zBinValues[i];

    // Check if track density already exists at current z position
    auto findIter =
        std::find(densityMap.zBin.begin(), densityMap.zBin.end(), z);

    if (findIter != densityMap.zBin.end()) {
      // Z bin already exists
      densityMap.density[std::distance(densityMap.zBin.begin(), findIter)] +=
          trackDensityMap.density[i];
    } else {
      // Create new z bin
      auto it =
          std::upper_bound(densityMap.zBin.begin(), densityMap.zBin.end(), z);
      densityMap.density.insert(densityMap.density.begin() +
                                    std::distance(densityMap.zBin.begin(), it),
                                trackDensityMap.density[i]);
      densityMap.zBin.insert(it, z);
    }
  }

  return trackDensityMap;
}

template <int trkGridSize>
void AdaptiveGridTrackDensity<trkGridSize>::subtractTrack(
    const TrackDensityMap& trackDensityMap,
    MainDensityMap& mainDensityMap) const {
  // Find position of current z bin in mainGridZValues
  auto findIter = std::find(mainDensityMap.zBin.begin(),
                            mainDensityMap.zBin.end(), trackDensityMap.zBin);
  // Calculate corresponding index in mainGridDensity
  std::int32_t densityIdx =
      std::distance(mainDensityMap.zBin.begin(), findIter);

  // Go over trkGrid and remove it from mainDensityGrid
  std::int32_t startEnd = static_cast<std::int32_t>((trkGridSize - 1) / 2);
  for (std::int32_t i = 0; i < trkGridSize; ++i) {
    mainDensityMap
        .density[static_cast<std::int32_t>(densityIdx + (i - startEnd))] -=
        trackDensityMap.density[i];
  }
}

template <int trkGridSize>
typename AdaptiveGridTrackDensity<trkGridSize>::TrackDensityMap
AdaptiveGridTrackDensity<trkGridSize>::createTrackGrid(
    std::int32_t centralZBin, const Vector2& impactParams,
    const SquareMatrix2& cov) const {
  TrackDensityMap trackGrid;

  trackGrid.zBin = centralZBin;

  for (std::int32_t i = 0; i < trkGridSize; ++i) {
    double z =
        (i - static_cast<double>(trkGridSize) / 2 + 0.5f) * m_cfg.binSize;
    // Bin coordinates in the d-z plane
    Vector2 binCoords(0., z);
    // Transformation to coordinate system with origin at the track center
    binCoords -= impactParams;
    trackGrid.density[i] = multivariateGaussian<2>(binCoords, cov);
  }

  return trackGrid;
}

template <int trkGridSize>
Result<double> AdaptiveGridTrackDensity<trkGridSize>::estimateSeedWidth(
    const MainDensityMap& mainDensityMap, double maxZ) const {
  if (mainDensityMap.empty()) {
    return VertexingError::EmptyInput;
  }
  // Get z bin of max density z value
  int sign = (maxZ > 0) ? +1 : -1;
  int zMaxGridBin = int(maxZ / m_cfg.binSize - sign * 0.5f);

  // Find location in mainGridZValues
  auto findIter = std::find(mainDensityMap.zBin.begin(),
                            mainDensityMap.zBin.end(), zMaxGridBin);
  int zBin = std::distance(mainDensityMap.zBin.begin(), findIter);

  const double maxValue = mainDensityMap.density[zBin];
  double gridValue = mainDensityMap.density[zBin];

  // Find right half-maximum bin
  int rhmBin = zBin;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if ((zMaxGridBin + (rhmBin - zBin)) != mainDensityMap.zBin[rhmBin]) {
      break;
    }
    rhmBin += 1;
    if (rhmBin == int(mainDensityMap.density.size())) {
      break;
    }
    gridValue = mainDensityMap.density[rhmBin];
  }

  // Use linear approximation to find better z value for FWHM between bins
  double deltaZ1 = (maxValue / 2 - mainDensityMap.density[rhmBin - 1]) *
                   (m_cfg.binSize / (mainDensityMap.density[rhmBin - 1] -
                                     mainDensityMap.density[rhmBin]));
  // Find left half-maximum bin
  int lhmBin = zBin;
  gridValue = mainDensityMap.density[zBin];
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if ((zMaxGridBin + (lhmBin - zBin)) != mainDensityMap.zBin[lhmBin]) {
      break;
    }
    lhmBin -= 1;
    if (lhmBin < 0) {
      break;
    }
    gridValue = mainDensityMap.density[lhmBin];
  }

  // Use linear approximation to find better z value for FWHM between bins
  double deltaZ2 = (maxValue / 2 - mainDensityMap.density[lhmBin + 1]) *
                   (m_cfg.binSize / (mainDensityMap.density[rhmBin + 1] -
                                     mainDensityMap.density[rhmBin]));

  // Approximate FWHM
  double fwhm =
      rhmBin * m_cfg.binSize - deltaZ1 - lhmBin * m_cfg.binSize - deltaZ2;

  // FWHM = 2.355 * sigma
  double width = fwhm / 2.355f;

  return std::isnormal(width) ? width : 0.0f;
}

template <int trkGridSize>
std::int32_t AdaptiveGridTrackDensity<trkGridSize>::highestDensitySumBin(
    MainDensityMap& mainDensityMap) const {
  // Checks the first (up to) 3 density maxima, if they are close, checks which
  // one has the highest surrounding density sum (the two neighboring bins)

  // The global maximum
  std::size_t zGridPos =
      std::distance(mainDensityMap.density.begin(),
                    std::max_element(mainDensityMap.density.begin(),
                                     mainDensityMap.density.end()));

  std::size_t zFirstMax = zGridPos;
  double firstDensity = mainDensityMap.density[zFirstMax];
  double firstSum = getDensitySum(mainDensityMap, zFirstMax);

  // Get the second highest maximum
  mainDensityMap.density[zFirstMax] = 0;
  zGridPos = std::distance(mainDensityMap.density.begin(),
                           std::max_element(mainDensityMap.density.begin(),
                                            mainDensityMap.density.end()));
  std::size_t zSecondMax = zGridPos;
  double secondDensity = mainDensityMap.density[zSecondMax];
  double secondSum = 0;
  if (firstDensity - secondDensity <
      firstDensity * m_cfg.maxRelativeDensityDev) {
    secondSum = getDensitySum(mainDensityMap, zSecondMax);
  }

  // Get the third highest maximum
  mainDensityMap.density[zSecondMax] = 0;
  zGridPos = std::distance(mainDensityMap.density.begin(),
                           std::max_element(mainDensityMap.density.begin(),
                                            mainDensityMap.density.end()));
  std::size_t zThirdMax = zGridPos;
  double thirdDensity = mainDensityMap.density[zThirdMax];
  double thirdSum = 0;
  if (firstDensity - thirdDensity <
      firstDensity * m_cfg.maxRelativeDensityDev) {
    thirdSum = getDensitySum(mainDensityMap, zThirdMax);
  }

  // Revert back to original values
  mainDensityMap.density[zFirstMax] = firstDensity;
  mainDensityMap.density[zSecondMax] = secondDensity;

  // Return the z-bin position of the highest density sum
  if (secondSum > firstSum && secondSum > thirdSum) {
    return zSecondMax;
  }
  if (thirdSum > secondSum && thirdSum > firstSum) {
    return zThirdMax;
  }
  return zFirstMax;
}

template <int trkGridSize>
double AdaptiveGridTrackDensity<trkGridSize>::getDensitySum(
    const MainDensityMap& mainDensityMap, std::size_t zBinIndex) const {
  double sum = mainDensityMap.density[zBinIndex];
  // Sum up only the density contributions from the neighboring bins if they are
  // still within bounds
  if (zBinIndex > 0) {
    // Check if we are still operating on continuous z values
    if (mainDensityMap.zBin[zBinIndex] - mainDensityMap.zBin[zBinIndex - 1] ==
        1) {
      sum += mainDensityMap.density[zBinIndex - 1];
    }
  }
  if (zBinIndex < mainDensityMap.size() - 1) {
    // Check if we are still operating on continuous z values
    if (mainDensityMap.zBin[zBinIndex + 1] - mainDensityMap.zBin[zBinIndex] ==
        1) {
      sum += mainDensityMap.density[zBinIndex + 1];
    }
  }
  return sum;
}

}  // namespace Acts
