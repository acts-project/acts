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
float Acts::AdaptiveGridTrackDensity<
    spatialTrkGridSize, temporalTrkGridSize>::getBinCenter(int bin,
                                                           float binExtent) {
  return bin * binExtent;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
int Acts::AdaptiveGridTrackDensity<
    spatialTrkGridSize, temporalTrkGridSize>::getBin(float value,
                                                     float binExtent) {
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
Acts::Result<typename Acts::AdaptiveGridTrackDensity<
    spatialTrkGridSize, temporalTrkGridSize>::ZTPosition>
Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    getMaxZTPosition(DensityMap& densityMap) const {
  if (densityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  Bin bin;
  if (!m_cfg.useHighestSumZPosition) {
    bin = highestDensityEntry(densityMap)->first;
  } else {
    // Get z position with highest density sum
    // of surrounding bins
    bin = highestDensitySumBin(densityMap);
  }

  // Derive corresponding z value
  float maxZ = getBinCenter(bin.first, m_cfg.spatialBinExtent);

  ZTPosition maxValues = std::make_pair(maxZ, 0.);

  // Get t value of the maximum if we do time vertex seeding
  if constexpr (temporalTrkGridSize > 1) {
    float maxT = getBinCenter(bin.second, m_cfg.temporalBinExtent);
    maxValues.second = maxT;
  }

  return maxValues;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
Acts::Result<typename Acts::AdaptiveGridTrackDensity<
    spatialTrkGridSize, temporalTrkGridSize>::ZTPositionAndWidth>
Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    getMaxZTPositionAndWidth(DensityMap& densityMap) const {
  // Get z value where the density is the highest
  auto maxZTRes = getMaxZTPosition(densityMap);
  if (not maxZTRes.ok()) {
    return maxZTRes.error();
  }
  ZTPosition maxZT = *maxZTRes;

  // Get seed width estimate
  auto widthRes = estimateSeedWidth(densityMap, maxZT);
  if (not widthRes.ok()) {
    return widthRes.error();
  }
  float width = *widthRes;
  ZTPositionAndWidth maxZTAndWidth{maxZT, width};
  return maxZTAndWidth;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
typename Acts::AdaptiveGridTrackDensity<spatialTrkGridSize,
                                        temporalTrkGridSize>::DensityMap
Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    addTrack(const Acts::BoundTrackParameters& trk,
             DensityMap& mainDensityMap) const {
  ActsVector<3> impactParams = trk.impactParameters();
  ActsSquareMatrix<3> cov = trk.impactParameterCovariance().value();

  // Calculate bin in d direction
  int centralDBin = getBin(impactParams(0), m_cfg.spatialBinExtent);
  // Check if current track affects grid density
  if (std::abs(centralDBin) > (spatialTrkGridSize - 1) / 2.) {
    DensityMap emptyTrackDensityMap;
    return emptyTrackDensityMap;
  }
  // Calculate bin in z direction
  int centralZBin = getBin(impactParams(1), m_cfg.spatialBinExtent);

  // If we don't do time vertex seeding, the time index is set to 0
  Bin centralBin = std::make_pair(centralZBin, 0.);

  // Calculate bin in t direction if we do time vertex seeding
  if constexpr (temporalTrkGridSize > 1) {
    int centralTBin = getBin(impactParams(2), m_cfg.temporalBinExtent);
    centralBin.second = centralTBin;
  }

  DensityMap trackDensityMap = createTrackGrid(impactParams, centralBin, cov);

  for (const auto& densityEntry : trackDensityMap) {
    Bin bin = densityEntry.first;
    float trackDensity = densityEntry.second;
    // Check if z bin is already part of the main grid
    if (mainDensityMap.count(bin) == 1) {
      mainDensityMap.at(bin) += trackDensity;
    } else {
      mainDensityMap[bin] = trackDensity;
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
    createTrackGrid(const Acts::Vector3& impactParams, const Bin& centralBin,
                    const Acts::SquareMatrix3& cov) const {
  DensityMap trackDensityMap;

  int halfSpatialTrkGridSize = (spatialTrkGridSize - 1) / 2;
  int firstZBin = centralBin.first - halfSpatialTrkGridSize;

  // If we don't do time vertex seeding, firstTBin will be 0.
  int halfTemporalTrkGridSize = (temporalTrkGridSize - 1) / 2;
  int firstTBin = centralBin.second - halfTemporalTrkGridSize;

  // Loop over bins
  for (int i = 0; i < temporalTrkGridSize; i++) {
    int tBin = firstTBin + i;
    // If we don't do vertex time seeding, we set the time to 0 since it will be
    // discarded in the for loop below anyways
    float t = 0;
    if constexpr (temporalTrkGridSize > 1) {
      t = getBinCenter(tBin, m_cfg.temporalBinExtent);
    }
    for (int j = 0; j < spatialTrkGridSize; j++) {
      int zBin = firstZBin + j;
      float z = getBinCenter(zBin, m_cfg.spatialBinExtent);
      // Bin coordinates in the d-z-t plane
      Acts::Vector3 binCoords(0., z, t);
      // Transformation to coordinate system with origin at the track center
      binCoords -= impactParams;
      Bin bin = std::make_pair(zBin, tBin);
      if constexpr (temporalTrkGridSize == 1) {
        trackDensityMap[bin] = multivariateGaussian<2>(
            binCoords.head<2>(), cov.topLeftCorner<2, 2>());
      } else {
        trackDensityMap[bin] = multivariateGaussian<3>(binCoords, cov);
      }
    }
  }
  return trackDensityMap;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
Acts::Result<float> Acts::AdaptiveGridTrackDensity<
    spatialTrkGridSize,
    temporalTrkGridSize>::estimateSeedWidth(const DensityMap& densityMap,
                                            const ZTPosition& maxZT) const {
  if (densityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  // Get z bin of max density
  int zMaxBin = getBin(maxZT.first, m_cfg.spatialBinExtent);
  int tMaxBin = 0;
  // Fill the time bin with a non-zero value if we do time vertex seeding
  if constexpr (temporalTrkGridSize > 1) {
    tMaxBin = getBin(maxZT.second, m_cfg.temporalBinExtent);
  }
  const float maxValue = densityMap.at(std::make_pair(zMaxBin, tMaxBin));

  int rhmBin = zMaxBin;
  float gridValue = maxValue;
  // Boolean indicating whether we find a filled bin that has a densityValue <=
  // maxValue/2
  bool binFilled = true;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (densityMap.count(std::make_pair(rhmBin + 1, tMaxBin)) == 0) {
      binFilled = false;
      break;
    }
    rhmBin += 1;
    gridValue = densityMap.at(std::make_pair(rhmBin, tMaxBin));
  }

  // Use linear approximation to find better z value for FWHM between bins
  float rightDensity = 0;
  if (binFilled) {
    rightDensity = densityMap.at(std::make_pair(rhmBin, tMaxBin));
  }
  float leftDensity = densityMap.at(std::make_pair(rhmBin - 1, tMaxBin));
  float deltaZ1 = m_cfg.spatialBinExtent * (maxValue / 2 - leftDensity) /
                  (rightDensity - leftDensity);

  int lhmBin = zMaxBin;
  gridValue = maxValue;
  binFilled = true;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (densityMap.count(std::make_pair(lhmBin - 1, tMaxBin)) == 0) {
      binFilled = false;
      break;
    }
    lhmBin -= 1;
    gridValue = densityMap.at(std::make_pair(lhmBin, tMaxBin));
  }

  // Use linear approximation to find better z value for FWHM between bins
  rightDensity = densityMap.at(std::make_pair(lhmBin + 1, tMaxBin));
  if (binFilled) {
    leftDensity = densityMap.at(std::make_pair(lhmBin, tMaxBin));
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
                         const Acts::ActsSquareMatrix<nDim>& cov) {
  float expo = -0.5 * args.transpose().dot(cov.inverse() * args);
  float gaussianDensity = safeExp(expo) / std::sqrt(cov.determinant());
  return gaussianDensity;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
typename Acts::AdaptiveGridTrackDensity<spatialTrkGridSize,
                                        temporalTrkGridSize>::Bin
Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    highestDensitySumBin(DensityMap& densityMap) const {
  // The global maximum
  auto firstMax = highestDensityEntry(densityMap);
  Bin binFirstMax = firstMax->first;
  float valueFirstMax = firstMax->second;
  float firstSum = getDensitySum(densityMap, binFirstMax);
  // Smaller maxima must have a density of at least:
  // valueFirstMax - densityDeviation
  float densityDeviation = valueFirstMax * m_cfg.maxRelativeDensityDev;

  // Get the second highest maximum
  densityMap.at(binFirstMax) = 0;
  auto secondMax = highestDensityEntry(densityMap);
  Bin binSecondMax = secondMax->first;
  float valueSecondMax = secondMax->second;
  float secondSum = 0;
  if (valueFirstMax - valueSecondMax < densityDeviation) {
    secondSum = getDensitySum(densityMap, binSecondMax);
  } else {
    // If the second maximum is not sufficiently large the third maximum won't
    // be either
    densityMap.at(binFirstMax) = valueFirstMax;
    return binFirstMax;
  }

  // Get the third highest maximum
  densityMap.at(binSecondMax) = 0;
  auto thirdMax = highestDensityEntry(densityMap);
  Bin binThirdMax = thirdMax->first;
  float valueThirdMax = thirdMax->second;
  float thirdSum = 0;
  if (valueFirstMax - valueThirdMax < densityDeviation) {
    thirdSum = getDensitySum(densityMap, binThirdMax);
  }

  // Revert back to original values
  densityMap.at(binFirstMax) = valueFirstMax;
  densityMap.at(binSecondMax) = valueSecondMax;

  // Return the z bin position of the highest density sum
  if (secondSum > firstSum && secondSum > thirdSum) {
    return binSecondMax;
  }
  if (thirdSum > secondSum && thirdSum > firstSum) {
    return binThirdMax;
  }
  return binFirstMax;
}

template <int spatialTrkGridSize, int temporalTrkGridSize>
float Acts::AdaptiveGridTrackDensity<spatialTrkGridSize, temporalTrkGridSize>::
    getDensitySum(const DensityMap& densityMap, const Bin& bin) const {
  // Add density from the bin.
  float sum = densityMap.at(bin);
  // Check if neighboring bins are part of the densityMap and add them (if they
  // are not part of the map, we assume them to be 0).
  // Note that each key in a map is unique; the .count() function therefore
  // returns either 0 or 1.
  Bin binShifted = bin;
  // Add density from the neighboring bin in -z direction.
  binShifted.first -= 1;
  if (densityMap.count(binShifted) == 1) {
    sum += densityMap.at(binShifted);
  }

  // Add density from the neighboring bin in +z direction.
  binShifted.first += 2;
  if (densityMap.count(binShifted) == 1) {
    sum += densityMap.at(binShifted);
  }
  return sum;
}
