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
double AdaptiveGridTrackDensity<trkGridSize>::getBinCenter(
    std::int32_t bin) const {
  return bin * m_cfg.binSize;
}

template <int trkGridSize>
std::int32_t AdaptiveGridTrackDensity<trkGridSize>::getBin(double value) const {
  return static_cast<std::int32_t>(std::floor(value / m_cfg.binSize - 0.5) + 1);
}

template <int trkGridSize>
Result<double> AdaptiveGridTrackDensity<trkGridSize>::getMaxZPosition(
    MainDensityMap& mainDensityMap) const {
  if (mainDensityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  std::size_t zBinIndex = 0;
  if (!m_cfg.useHighestSumZPosition) {
    zBinIndex = std::distance(mainDensityMap.density.begin(),
                              std::max_element(mainDensityMap.density.begin(),
                                               mainDensityMap.density.end()));
  } else {
    // Get z position with highest density sum of surrounding bins
    zBinIndex = highestDensitySumBinIndex(mainDensityMap);
  }

  std::int32_t zBin = mainDensityMap.zBin[zBinIndex];
  return getBinCenter(zBin);
}

template <int trkGridSize>
Result<typename AdaptiveGridTrackDensity<trkGridSize>::ZPositionAndWidth>
AdaptiveGridTrackDensity<trkGridSize>::getMaxZPositionAndWidth(
    MainDensityMap& mainDensityMap) const {
  // Get z maximum value
  auto maxZRes = getMaxZPosition(mainDensityMap);
  if (!maxZRes.ok()) {
    return maxZRes.error();
  }
  double maxZ = *maxZRes;

  // Get seed width estimate
  auto widthRes = estimateSeedWidth(mainDensityMap, maxZ);
  if (!widthRes.ok()) {
    return widthRes.error();
  }
  double width = *widthRes;

  return ZPositionAndWidth{maxZ, width};
}

template <int trkGridSize>
typename AdaptiveGridTrackDensity<trkGridSize>::TrackDensityMap
AdaptiveGridTrackDensity<trkGridSize>::addTrack(
    const BoundTrackParameters& trk, MainDensityMap& mainDensityMap) const {
  ActsVector<2> impactParams = trk.spatialImpactParameters();
  ActsSquareMatrix<2> cov = trk.spatialImpactParameterCovariance().value();

  double d0 = trk.parameters()[0];
  double z0 = trk.parameters()[1];

  // Calculate bin in d direction
  std::int32_t centralDBin = getBin(d0);

  std::int32_t halfTrkGridSize = (trkGridSize - 1) / 2;

  // Check if current track affects grid density
  if (std::abs(centralDBin) > halfTrkGridSize) {
    return {};
  }

  // Calculate bin in z
  std::int32_t centralZBin = getBin(z0);

  TrackDensityMap trackDensityMap =
      createTrackGrid(centralZBin, impactParams, cov);

  std::int32_t firstZBin = centralZBin - halfTrkGridSize;

  for (std::int32_t i = 0; i < trkGridSize; ++i) {
    std::int32_t zBin = firstZBin + i;

    // Check if track density already exists at current z position
    if (auto findIter = std::find(mainDensityMap.zBin.begin(),
                                  mainDensityMap.zBin.end(), zBin);
        findIter != mainDensityMap.zBin.end()) {
      // z bin already exists
      std::size_t zBinIndex =
          std::distance(mainDensityMap.zBin.begin(), findIter);
      mainDensityMap.density[zBinIndex] += trackDensityMap.density[i];
    } else {
      // Create new z bin
      auto it = std::upper_bound(mainDensityMap.zBin.begin(),
                                 mainDensityMap.zBin.end(), zBin);
      mainDensityMap.density.insert(
          mainDensityMap.density.begin() +
              std::distance(mainDensityMap.zBin.begin(), it),
          trackDensityMap.density[i]);
      mainDensityMap.zBin.insert(it, zBin);
    }
  }

  return trackDensityMap;
}

template <int trkGridSize>
void AdaptiveGridTrackDensity<trkGridSize>::subtractTrack(
    const TrackDensityMap& trackDensityMap,
    MainDensityMap& mainDensityMap) const {
  // Calculate corresponding index in mainDensityMap
  std::size_t centralZIndex = std::distance(
      mainDensityMap.zBin.begin(),
      std::find(mainDensityMap.zBin.begin(), mainDensityMap.zBin.end(),
                trackDensityMap.centralZBin));
  std::int32_t halfTrkGridSize = (trkGridSize - 1) / 2;
  std::int32_t firstZBin = centralZIndex - halfTrkGridSize;

  // Go over trkGrid and remove it from mainDensityGrid
  for (std::int32_t i = 0; i < trkGridSize; ++i) {
    mainDensityMap.density[firstZBin + i] -= trackDensityMap.density[i];
  }
}

template <int trkGridSize>
typename AdaptiveGridTrackDensity<trkGridSize>::TrackDensityMap
AdaptiveGridTrackDensity<trkGridSize>::createTrackGrid(
    std::int32_t centralZBin, const Vector2& impactParams,
    const SquareMatrix2& cov) const {
  TrackDensityMap trackGrid;
  trackGrid.centralZBin = centralZBin;

  std::int32_t halfTrkGridSize = (trkGridSize - 1) / 2;
  std::int32_t firstZBin = centralZBin - halfTrkGridSize;
  for (std::int32_t i = 0; i < trkGridSize; ++i) {
    std::int32_t zBin = firstZBin + i;
    double z = getBinCenter(zBin);
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
  std::int32_t zMaxBin = getBin(maxZ);

  // Find location in mainDensityMap
  std::size_t zMaxIndex =
      std::distance(mainDensityMap.zBin.begin(),
                    std::find(mainDensityMap.zBin.begin(),
                              mainDensityMap.zBin.end(), zMaxBin));

  const double maxValue = mainDensityMap.density[zMaxIndex];

  // Find right half-maximum bin
  std::size_t rhmIndex = zMaxIndex;
  double gridValue = mainDensityMap.density[zMaxIndex];
  // Boolean indicating whether we find a filled bin that has a densityValue <=
  // maxValue/2
  bool binFilled = true;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if ((zMaxBin + static_cast<std::int32_t>(rhmIndex - zMaxIndex)) !=
        mainDensityMap.zBin[rhmIndex]) {
      binFilled = false;
      break;
    }
    ++rhmIndex;
    if (rhmIndex == mainDensityMap.density.size()) {
      break;
    }
    gridValue = mainDensityMap.density[rhmIndex];
  }

  // Use linear approximation to find better z value for FWHM between bins
  double rightDensity = 0;
  if (binFilled) {
    rightDensity = mainDensityMap.density[rhmIndex];
  }
  double leftDensity = mainDensityMap.density[rhmIndex - 1];
  double deltaZ1 = m_cfg.binSize * (maxValue / 2 - leftDensity) /
                   (rightDensity - leftDensity);

  // Find left half-maximum bin
  std::size_t lhmIndex = zMaxIndex;
  gridValue = maxValue;
  binFilled = true;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if ((zMaxBin + static_cast<std::int32_t>(lhmIndex - zMaxIndex)) !=
        mainDensityMap.zBin[lhmIndex]) {
      binFilled = false;
      break;
    }
    if (lhmIndex == 0) {
      break;
    }
    --lhmIndex;
    gridValue = mainDensityMap.density[lhmIndex];
  }

  // Use linear approximation to find better z value for FWHM between bins
  rightDensity = mainDensityMap.density[lhmIndex + 1];
  if (binFilled) {
    leftDensity = mainDensityMap.density[lhmIndex];
  } else {
    leftDensity = 0;
  }
  double deltaZ2 = m_cfg.binSize * (rightDensity - maxValue / 2) /
                   (rightDensity - leftDensity);

  // Approximate FWHM
  double fwhm = static_cast<double>(rhmIndex) * m_cfg.binSize - deltaZ1 -
                static_cast<double>(lhmIndex) * m_cfg.binSize - deltaZ2;

  // FWHM = 2.355 * sigma
  double width = fwhm / 2.355f;

  return std::isnormal(width) ? width : 0.0f;
}

template <int trkGridSize>
std::size_t AdaptiveGridTrackDensity<trkGridSize>::highestDensitySumBinIndex(
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
