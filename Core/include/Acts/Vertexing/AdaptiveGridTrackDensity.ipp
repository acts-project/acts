// This file is part of the Acts project.
//
// Copyright (C) 2020-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>

template <int spatialTrkGridSize, int temporalTrkGridSize>
float Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    getBinCenter(int bin, float binExtent) const {
  return bin * binExtent;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
int Acts::AdaptiveGridTrackDensity<
    spatialTrkGridSize, temporalTrkGridSize>::getBin(float value,
                                                     float binExtent) const {
  return static_cast<int>(std::floor(value / binExtent - 0.5) + 1);
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
typename Acts::AdaptiveGridTrackDensity<
    spatialTrkGridSize, temporalTrkGridSize>::DensityMap::const_iterator
Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    highestDensityEntry(const DensityMap& densityMap) const {
  auto maxEntry = std::max_element(
      std::begin(densityMap), std::end(densityMap),
      [](const auto& densityEntry1, const auto& densityEntry2) {
        return densityEntry1.second < densityEntry2.second;
      });
  return maxEntry;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
Acts::Result<float> Acts::AdaptiveGridTrackDensity<
    spatialTrkGridSize,
    temporalTrkGridSize>::getMaxZPosition(DensityMap& densityMap) const {
  if (densityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  int zBin = -1;
  if (!m_cfg.useHighestSumZPosition) {
    zBin = highestDensityEntry(densityMap)->first;
  } else {
    // Get z position with highest density sum
    // of surrounding bins
    zBin = highestDensitySumBin(densityMap);
  }

  // Derive corresponding z value
  return getBinCenter(zBin, m_cfg.spatialBinExtent);
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
Acts::Result<std::pair<float, float>>
Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    getMaxZPositionAndWidth(DensityMap& densityMap) const {
  // Get z value where the density is the highest
  auto maxZRes = getMaxZPosition(densityMap);
  if (not maxZRes.ok()) {
    return maxZRes.error();
  }
  float maxZ = *maxZRes;

  // Get seed width estimate
  auto widthRes = estimateSeedWidth(densityMap, maxZ);
  if (not widthRes.ok()) {
    return widthRes.error();
  }
  float width = *widthRes;
  std::pair<float, float> returnPair{maxZ, width};
  return returnPair;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
typename Acts::AdaptiveGridTrackDensity<spatialTrkGridSize,
                                        temporalTrkGridSize>::DensityMap
Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    addTrack(const Acts::BoundTrackParameters& trk,
             DensityMap& mainDensityMap) const {
  SquareMatrix2 cov = trk.covariance().value().block<2, 2>(0, 0);
  ActsVector<2> impactParams = trk.parameters().head<2>();

  // Calculate bin in d direction
  int centralDBin = getBin(impactParams(0), m_cfg.spatialBinExtent);
  // Check if current track affects grid density
  if (std::abs(centralDBin) > (spatialTrkGridSize - 1) / 2.) {
    DensityMap emptyTrackDensityMap;
    return emptyTrackDensityMap;
  }
  // Calculate bin in z direction
  int centralZBin = getBin(impactParams(1), m_cfg.spatialBinExtent);

  DensityMap trackDensityMap = createTrackGrid(impactParams, centralZBin, cov);

  for (const auto& densityEntry : trackDensityMap) {
    int zBin = densityEntry.first;
    float trackDensity = densityEntry.second;
    // Check if z bin is already part of the main grid
    if (mainDensityMap.count(zBin) == 1) {
      mainDensityMap.at(zBin) += trackDensity;
    } else {
      mainDensityMap[zBin] = trackDensity;
    }
  }

  return trackDensityMap;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
void Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    subtractTrack(const DensityMap& trackDensityMap,
                  DensityMap& mainDensityMap) const {
  for (auto it = trackDensityMap.begin(); it != trackDensityMap.end(); it++) {
    mainDensityMap.at(it->first) -= it->second;
  }
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
typename Acts::AdaptiveGridTrackDensity<spatialTrkGridSize,
                                        temporalTrkGridSize>::DensityMap
Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    createTrackGrid(const Acts::ActsVector<2>& impactParams, int centralZBin,
                    const Acts::SquareMatrix2& cov) const {
  DensityMap trackDensityMap;

  int halfTrkGridSize = (spatialTrkGridSize - 1) / 2;
  int firstZBin = centralZBin - halfTrkGridSize;
  // Loop over bins
  for (int j = 0; j < spatialTrkGridSize; j++) {
    int zBin = firstZBin + j;
    float z = getBinCenter(zBin, m_cfg.spatialBinExtent);
    // Transverse and logitudinal coordinate of the bin
    Acts::Vector2 binCoords(0., z);
    // Transformation to coordinate system with origin at the track center
    binCoords -= impactParams;
    trackDensityMap[zBin] = multivariateGaussian<2>(binCoords, cov);
  }
  return trackDensityMap;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
Acts::Result<float> Acts::AdaptiveGridTrackDensity<
    spatialTrkGridSize,
    temporalTrkGridSize>::estimateSeedWidth(const DensityMap& densityMap,
                                            float maxZ) const {
  if (densityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  // Get z bin of max density and max density value
  int zMaxBin = getBin(maxZ, m_cfg.spatialBinExtent);
  const float maxValue = densityMap.at(zMaxBin);

  int rhmBin = zMaxBin;
  float gridValue = maxValue;
  // Boolean indicating whether we find a filled bin that has a densityValue <=
  // maxValue/2
  bool binFilled = true;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (densityMap.count(rhmBin + 1) == 0) {
      binFilled = false;
      break;
    }
    rhmBin += 1;
    gridValue = densityMap.at(rhmBin);
  }

  // Use linear approximation to find better z value for FWHM between bins
  float rightDensity = 0;
  if (binFilled) {
    rightDensity = densityMap.at(rhmBin);
  }
  float leftDensity = densityMap.at(rhmBin - 1);
  float deltaZ1 = m_cfg.spatialBinExtent * (maxValue / 2 - leftDensity) /
                  (rightDensity - leftDensity);

  int lhmBin = zMaxBin;
  gridValue = maxValue;
  binFilled = true;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (densityMap.count(lhmBin - 1) == 0) {
      binFilled = false;
      break;
    }
    lhmBin -= 1;
    gridValue = densityMap.at(lhmBin);
  }

  // Use linear approximation to find better z value for FWHM between bins
  rightDensity = densityMap.at(lhmBin + 1);
  if (binFilled) {
    leftDensity = densityMap.at(lhmBin);
  } else {
    leftDensity = 0;
  }
  float deltaZ2 = m_cfg.spatialBinExtent * (rightDensity - maxValue / 2) /
                  (rightDensity - leftDensity);

  float fwhm = (rhmBin - lhmBin) * m_cfg.spatialBinExtent - deltaZ1 - deltaZ2;

  // FWHM = 2.355 * sigma
  float width = fwhm / 2.355f;

  return std::isnormal(width) ? width : 0.0f;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
template <unsigned int nDim>
float Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    multivariateGaussian(const Acts::ActsVector<nDim>& args,
                         const Acts::ActsSquareMatrix<nDim>& cov) const {
  ActsScalar expo = -0.5 * args.transpose().dot(cov.inverse() * args);
  // We perform the calculation in the ActsScalar precision but we return a
  // float
  ActsScalar gaussianDensity = safeExp(expo) / std::sqrt(cov.determinant());
  return static_cast<float>(gaussianDensity);
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
int Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    highestDensitySumBin(DensityMap& densityMap) const {
  // Checks the first (up to) 3 density maxima, if they are close, checks which
  // one has the highest surrounding density sum (the two neighboring bins)

  // The global maximum
  auto firstMax = highestDensityEntry(densityMap);
  int zBinFirstMax = firstMax->first;
  float valueFirstMax = firstMax->second;
  float firstSum = getDensitySum(densityMap, zBinFirstMax);
  float densityDeviation = valueFirstMax * m_cfg.maxRelativeDensityDev;

  // Get the second highest maximum
  densityMap.at(zBinFirstMax) = 0;
  auto secondMax = highestDensityEntry(densityMap);
  int zBinSecondMax = secondMax->first;
  float valueSecondMax = secondMax->second;
  float secondSum = 0;
  if (valueFirstMax - valueSecondMax < densityDeviation) {
    secondSum = getDensitySum(densityMap, zBinSecondMax);
  }

  // Get the third highest maximum
  densityMap.at(zBinSecondMax) = 0;
  auto thirdMax = highestDensityEntry(densityMap);
  int zBinThirdMax = thirdMax->first;
  float valueThirdMax = thirdMax->second;
  float thirdSum = 0;
  if (valueFirstMax - valueThirdMax < densityDeviation) {
    thirdSum = getDensitySum(densityMap, zBinThirdMax);
  }

  // Revert back to original values
  densityMap.at(zBinFirstMax) = valueFirstMax;
  densityMap.at(zBinSecondMax) = valueSecondMax;

  // Return the z bin position of the highest density sum
  if (secondSum > firstSum && secondSum > thirdSum) {
    return zBinSecondMax;
  }
  if (thirdSum > secondSum && thirdSum > firstSum) {
    return zBinThirdMax;
  }
  return zBinFirstMax;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
float Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    getDensitySum(const DensityMap& densityMap, int zBin) const {
  float sum = densityMap.at(zBin);
  // Check if neighboring bins are part of the densityMap and add them (if they
  // are not part of the map, we assume them to be 0) Note that each key in a
  // map is unique; the .count() function therefore returns either 0 or 1
  if (densityMap.count(zBin - 1) == 1) {
    sum += densityMap.at(zBin - 1);
  }
  if (densityMap.count(zBin + 1) == 1) {
    sum += densityMap.at(zBin + 1);
  }
  return sum;
}
