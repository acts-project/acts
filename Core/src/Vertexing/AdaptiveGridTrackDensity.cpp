// This file is part of the Acts project.
//
// Copyright (C) 2020-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "Acts/Vertexing/AdaptiveGridTrackDensity.hpp"

#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>

namespace Acts {

double AdaptiveGridTrackDensity::getBinCenter(int bin, double binExtent) {
  return bin * binExtent;
}

int AdaptiveGridTrackDensity::getBin(double value, double binExtent) {
  return static_cast<int>(std::floor(value / binExtent - 0.5) + 1);
}

int AdaptiveGridTrackDensity::getTrkGridSize(double sigma, double trkSigmas,
                                             double binExtent) {
  if (trkSigmas == 0) {
    return 1;
  }
  int size = static_cast<int>(std::ceil(2 * trkSigmas * sigma / binExtent));
  if (size % 2 == 0) {
    ++size;
  }
  return size;
}

int AdaptiveGridTrackDensity::getSpatialTrkGridSize(double sigma) const {
  return getTrkGridSize(sigma, m_cfg.spatialTrkSigmas, m_cfg.spatialBinExtent);
}

int AdaptiveGridTrackDensity::getTemporalTrkGridSize(double sigma) const {
  return getTrkGridSize(sigma, m_cfg.temporalTrkSigmas,
                        m_cfg.temporalBinExtent);
}

AdaptiveGridTrackDensity::AdaptiveGridTrackDensity(const Config& cfg)
    : m_cfg(cfg) {}

AdaptiveGridTrackDensity::DensityMap::const_iterator
AdaptiveGridTrackDensity::highestDensityEntry(
    const DensityMap& densityMap) const {
  auto maxEntry = std::max_element(
      std::begin(densityMap), std::end(densityMap),
      [](const auto& densityEntry1, const auto& densityEntry2) {
        return densityEntry1.second < densityEntry2.second;
      });
  return maxEntry;
}

Result<AdaptiveGridTrackDensity::ZTPosition>
AdaptiveGridTrackDensity::getMaxZTPosition(DensityMap& densityMap) const {
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

  // Derive corresponding z and t value
  double maxZ = getBinCenter(bin.first, m_cfg.spatialBinExtent);
  double maxT = getBinCenter(bin.second, m_cfg.temporalBinExtent);

  return std::make_pair(maxZ, maxT);
}

Result<AdaptiveGridTrackDensity::ZTPositionAndWidth>
AdaptiveGridTrackDensity::getMaxZTPositionAndWidth(
    DensityMap& densityMap) const {
  // Get z value where the density is the highest
  auto maxZTRes = getMaxZTPosition(densityMap);
  if (!maxZTRes.ok()) {
    return maxZTRes.error();
  }
  ZTPosition maxZT = *maxZTRes;

  // Get seed width estimate
  auto widthRes = estimateSeedWidth(densityMap, maxZT);
  if (!widthRes.ok()) {
    return widthRes.error();
  }
  double width = *widthRes;
  ZTPositionAndWidth maxZTAndWidth{maxZT, width};
  return maxZTAndWidth;
}

AdaptiveGridTrackDensity::DensityMap AdaptiveGridTrackDensity::addTrack(
    const BoundTrackParameters& trk, DensityMap& mainDensityMap) const {
  ActsVector<3> impactParams = trk.impactParameters();
  ActsSquareMatrix<3> cov = trk.impactParameterCovariance().value();

  int spatialTrkGridSize = getSpatialTrkGridSize(cov(1, 1));
  int temporalTrkGridSize = getTemporalTrkGridSize(cov(2, 2));

  // Calculate bin in d direction
  int centralDBin = getBin(impactParams(0), m_cfg.spatialBinExtent);
  // Check if current track affects grid density
  if (std::abs(centralDBin) > (spatialTrkGridSize - 1) / 2.) {
    // return an empty map
    return {};
  }

  // Calculate bin in z and t direction
  int centralZBin = getBin(impactParams(1), m_cfg.spatialBinExtent);
  int centralTBin = getBin(impactParams(2), m_cfg.temporalBinExtent);

  Bin centralBin = {centralZBin, centralTBin};

  DensityMap trackDensityMap = createTrackGrid(
      impactParams, centralBin, cov, spatialTrkGridSize, temporalTrkGridSize);

  for (const auto& [bin, trackDensity] : trackDensityMap) {
    mainDensityMap[bin] += trackDensity;
  }

  return trackDensityMap;
}

void AdaptiveGridTrackDensity::subtractTrack(const DensityMap& trackDensityMap,
                                             DensityMap& mainDensityMap) const {
  for (const auto& [bin, trackDensity] : trackDensityMap) {
    mainDensityMap[bin] -= trackDensity;
  }
}

AdaptiveGridTrackDensity::DensityMap AdaptiveGridTrackDensity::createTrackGrid(
    const Vector3& impactParams, const Bin& centralBin,
    const SquareMatrix3& cov, int spatialTrkGridSize,
    int temporalTrkGridSize) const {
  DensityMap trackDensityMap;

  int halfSpatialTrkGridSize = (spatialTrkGridSize - 1) / 2;
  int firstZBin = centralBin.first - halfSpatialTrkGridSize;

  // If we don't do time vertex seeding, firstTBin will be 0.
  int halfTemporalTrkGridSize = (temporalTrkGridSize - 1) / 2;
  int firstTBin = centralBin.second - halfTemporalTrkGridSize;

  // Loop over bins
  for (int i = 0; i < temporalTrkGridSize; i++) {
    int tBin = firstTBin + i;
    double t = getBinCenter(tBin, m_cfg.temporalBinExtent);
    if (t < m_cfg.temporalWindow.first || t > m_cfg.temporalWindow.second) {
      continue;
    }
    for (int j = 0; j < spatialTrkGridSize; j++) {
      int zBin = firstZBin + j;
      double z = getBinCenter(zBin, m_cfg.spatialBinExtent);
      if (z < m_cfg.spatialWindow.first || z > m_cfg.spatialWindow.second) {
        continue;
      }
      // Bin coordinates in the d-z-t plane
      Vector3 binCoords(0., z, t);
      // Transformation to coordinate system with origin at the track center
      binCoords -= impactParams;
      Bin bin = {zBin, tBin};
      trackDensityMap[bin] = multivariateGaussian<3>(binCoords, cov);
    }
  }

  return trackDensityMap;
}

Result<double> AdaptiveGridTrackDensity::estimateSeedWidth(
    const DensityMap& densityMap, const ZTPosition& maxZT) const {
  if (densityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  // Get z and t bin of max density
  int zMaxBin = getBin(maxZT.first, m_cfg.spatialBinExtent);
  int tMaxBin = getBin(maxZT.second, m_cfg.temporalBinExtent);
  const double maxValue = densityMap.at({zMaxBin, tMaxBin});

  int rhmBin = zMaxBin;
  double gridValue = maxValue;
  // Boolean indicating whether we find a filled bin that has a densityValue <=
  // maxValue/2
  bool binFilled = true;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (densityMap.count({rhmBin + 1, tMaxBin}) == 0) {
      binFilled = false;
      break;
    }
    rhmBin += 1;
    gridValue = densityMap.at({rhmBin, tMaxBin});
  }

  // Use linear approximation to find better z value for FWHM between bins
  double rightDensity = 0;
  if (binFilled) {
    rightDensity = densityMap.at({rhmBin, tMaxBin});
  }
  double leftDensity = densityMap.at({rhmBin - 1, tMaxBin});
  double deltaZ1 = m_cfg.spatialBinExtent * (maxValue / 2 - leftDensity) /
                   (rightDensity - leftDensity);

  int lhmBin = zMaxBin;
  gridValue = maxValue;
  binFilled = true;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (densityMap.count({lhmBin - 1, tMaxBin}) == 0) {
      binFilled = false;
      break;
    }
    lhmBin -= 1;
    gridValue = densityMap.at({lhmBin, tMaxBin});
  }

  // Use linear approximation to find better z value for FWHM between bins
  rightDensity = densityMap.at({lhmBin + 1, tMaxBin});
  if (binFilled) {
    leftDensity = densityMap.at({lhmBin, tMaxBin});
  } else {
    leftDensity = 0;
  }
  double deltaZ2 = m_cfg.spatialBinExtent * (rightDensity - maxValue / 2) /
                   (rightDensity - leftDensity);

  double fwhm = (rhmBin - lhmBin) * m_cfg.spatialBinExtent - deltaZ1 - deltaZ2;

  // FWHM = 2.355 * sigma
  double width = fwhm / 2.355f;

  return std::isnormal(width) ? width : 0.0f;
}

template <unsigned int nDim>
double AdaptiveGridTrackDensity::multivariateGaussian(
    const ActsVector<nDim>& args, const ActsSquareMatrix<nDim>& cov) {
  double expo = -0.5 * args.transpose().dot(cov.inverse() * args);
  double gaussianDensity = safeExp(expo) / std::sqrt(cov.determinant());
  return gaussianDensity;
}

AdaptiveGridTrackDensity::Bin AdaptiveGridTrackDensity::highestDensitySumBin(
    DensityMap& densityMap) const {
  // The global maximum
  auto firstMax = highestDensityEntry(densityMap);
  Bin binFirstMax = firstMax->first;
  double valueFirstMax = firstMax->second;
  double firstSum = getDensitySum(densityMap, binFirstMax);
  // Smaller maxima must have a density of at least:
  // valueFirstMax - densityDeviation
  double densityDeviation = valueFirstMax * m_cfg.maxRelativeDensityDev;

  // Get the second highest maximum
  densityMap.at(binFirstMax) = 0;
  auto secondMax = highestDensityEntry(densityMap);
  Bin binSecondMax = secondMax->first;
  double valueSecondMax = secondMax->second;
  double secondSum = 0;
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
  double valueThirdMax = thirdMax->second;
  double thirdSum = 0;
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

double AdaptiveGridTrackDensity::getDensitySum(const DensityMap& densityMap,
                                               const Bin& bin) const {
  auto valueOrDefault = [&densityMap](Bin b) -> double {
    auto it = densityMap.find(b);
    if (it == densityMap.end()) {
      return 0;
    }
    return it->second;
  };

  // Add density from the bin.
  // Check if neighboring bins are part of the densityMap and add them (if they
  // are not part of the map, we assume them to be 0).
  return valueOrDefault(bin) + valueOrDefault({bin.first, bin.second - 1}) +
         valueOrDefault({bin.first, bin.second + 1});
}

}  // namespace Acts
