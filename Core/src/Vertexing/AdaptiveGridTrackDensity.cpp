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

namespace {

template <unsigned int nDim>
double multivariateGaussianExponent(const ActsVector<nDim>& args,
                                    const ActsSquareMatrix<nDim>& cov) {
  return -0.5 * args.transpose().dot(cov.inverse() * args);
}

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
                            const ActsSquareMatrix<nDim>& cov,
                            double exponentOffset) {
  double exponent = multivariateGaussianExponent<nDim>(args, cov);
  double gaussianDensity =
      safeExp(exponent + exponentOffset) / std::sqrt(cov.determinant());
  return gaussianDensity;
}

double rescaleDensityMapFactor(double oldExponentOffset,
                               double newExponentOffset) {
  return safeExp(newExponentOffset - oldExponentOffset);
}

void rescaleDensityMap(AdaptiveGridTrackDensity::DensityMap& densityMap,
                       double newExponentOffset) {
  double rescaleFactor =
      rescaleDensityMapFactor(densityMap.exponentOffset, newExponentOffset);

  for (auto& [bin, density] : densityMap.scaledMap) {
    density *= rescaleFactor;
  }
  densityMap.exponentOffset = newExponentOffset;
}

}  // namespace

double AdaptiveGridTrackDensity::getBinCenter(int bin, double binExtent) {
  return bin * binExtent;
}

int AdaptiveGridTrackDensity::getBin(double value, double binExtent) {
  return static_cast<int>(std::floor(value / binExtent - 0.5) + 1);
}

int AdaptiveGridTrackDensity::getTrkGridSize(double sigma, double trkSigmas,
                                             double binExtent) {
  int size = static_cast<int>(std::ceil(2 * trkSigmas * sigma / binExtent));
  if (size % 2 == 0) {
    ++size;
  }
  return size;
}

int AdaptiveGridTrackDensity::getSpatialBin(double value) const {
  return getBin(value, m_cfg.spatialBinExtent);
}

int AdaptiveGridTrackDensity::getTemporalBin(double value) const {
  if (!m_cfg.useTime) {
    return 0;
  }
  return getBin(value, m_cfg.temporalBinExtent);
}

double AdaptiveGridTrackDensity::getSpatialBinCenter(int bin) const {
  return getBinCenter(bin, m_cfg.spatialBinExtent);
}

double AdaptiveGridTrackDensity::getTemporalBinCenter(int bin) const {
  if (!m_cfg.useTime) {
    return 0.;
  }
  return getBinCenter(bin, m_cfg.temporalBinExtent);
}

int AdaptiveGridTrackDensity::getSpatialTrkGridSize(double sigma) const {
  return getTrkGridSize(sigma, m_cfg.spatialTrkSigmas, m_cfg.spatialBinExtent);
}

int AdaptiveGridTrackDensity::getTemporalTrkGridSize(double sigma) const {
  if (!m_cfg.useTime) {
    return 1;
  }
  return getTrkGridSize(sigma, m_cfg.temporalTrkSigmas,
                        m_cfg.temporalBinExtent);
}

AdaptiveGridTrackDensity::AdaptiveGridTrackDensity(const Config& cfg)
    : m_cfg(cfg) {}

AdaptiveGridTrackDensity::SparseDensityMap::const_iterator
AdaptiveGridTrackDensity::highestDensityEntry(
    const DensityMap& densityMap) const {
  auto maxEntry = std::max_element(
      std::begin(densityMap.scaledMap), std::end(densityMap.scaledMap),
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
  double maxZ = getSpatialBinCenter(bin.first);
  double maxT = getTemporalBinCenter(bin.second);

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

  int spatialTrkGridSize = getSpatialTrkGridSize(std::sqrt(cov(1, 1)));
  int temporalTrkGridSize = getTemporalTrkGridSize(std::sqrt(cov(2, 2)));

  // Calculate bin in z and t direction
  int centralZBin = getSpatialBin(impactParams(1));
  int centralTBin = getTemporalBin(impactParams(2));

  Bin centralBin = {centralZBin, centralTBin};

  DensityMap trackDensityMap = createTrackGrid(
      impactParams, centralBin, cov, spatialTrkGridSize, temporalTrkGridSize);

  if (mainDensityMap.exponentOffset < trackDensityMap.exponentOffset) {
    double rescaleFactor = rescaleDensityMapFactor(
        trackDensityMap.exponentOffset, mainDensityMap.exponentOffset);

    for (const auto& [bin, density] : trackDensityMap.scaledMap) {
      mainDensityMap.scaled(bin) += density * rescaleFactor;
    }
  } else {
    rescaleDensityMap(mainDensityMap, trackDensityMap.exponentOffset);

    for (const auto& [bin, density] : trackDensityMap.scaledMap) {
      mainDensityMap.scaled(bin) += density;
    }
  }

  return trackDensityMap;
}

void AdaptiveGridTrackDensity::subtractTrack(const DensityMap& trackDensityMap,
                                             DensityMap& mainDensityMap) const {
  double rescaleFactor = rescaleDensityMapFactor(trackDensityMap.exponentOffset,
                                                 mainDensityMap.exponentOffset);

  for (const auto& [bin, density] : trackDensityMap.scaledMap) {
    mainDensityMap.scaled(bin) -= density * rescaleFactor;
  }
}

AdaptiveGridTrackDensity::DensityMap AdaptiveGridTrackDensity::createTrackGrid(
    const Vector3& impactParams, const Bin& centralBin,
    const SquareMatrix3& cov, int spatialTrkGridSize,
    int temporalTrkGridSize) const {
  DensityMap trackDensityMap;

  Vector3 maxParams(impactParams[0], 0, 0);
  double exponentOffset = 0;
  if (m_cfg.useTime) {
    exponentOffset -= multivariateGaussianExponent<3>(maxParams, cov);
  } else {
    exponentOffset -= multivariateGaussianExponent<2>(
        maxParams.head<2>(), cov.topLeftCorner<2, 2>());
  }
  trackDensityMap.exponentOffset = exponentOffset;

  int halfSpatialTrkGridSize = (spatialTrkGridSize - 1) / 2;
  int firstZBin = centralBin.first - halfSpatialTrkGridSize;

  // If we don't do time vertex seeding, firstTBin will be 0.
  int halfTemporalTrkGridSize = (temporalTrkGridSize - 1) / 2;
  int firstTBin = centralBin.second - halfTemporalTrkGridSize;

  // Loop over bins
  for (int i = 0; i < temporalTrkGridSize; i++) {
    int tBin = firstTBin + i;
    double t = getTemporalBinCenter(tBin);
    if (t < m_cfg.temporalWindow.first || t > m_cfg.temporalWindow.second) {
      continue;
    }
    for (int j = 0; j < spatialTrkGridSize; j++) {
      int zBin = firstZBin + j;
      double z = getSpatialBinCenter(zBin);
      if (z < m_cfg.spatialWindow.first || z > m_cfg.spatialWindow.second) {
        continue;
      }
      // Bin coordinates in the d-z-t plane
      Vector3 binCoords(0., z, t);
      // Transformation to coordinate system with origin at the track center
      binCoords -= impactParams;
      Bin bin = {zBin, tBin};
      double density = 0;
      if (m_cfg.useTime) {
        density = multivariateGaussian<3>(binCoords, cov, exponentOffset);
      } else {
        density = multivariateGaussian<2>(
            binCoords.head<2>(), cov.topLeftCorner<2, 2>(), exponentOffset);
      }
      if (density > 0) {
        trackDensityMap.scaled(bin) = density;
      }
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
  double maxValue = densityMap.scaled({zMaxBin, tMaxBin});

  int rhmBin = zMaxBin;
  double gridValue = maxValue;
  // Boolean indicating whether we find a filled bin that has a densityValue <=
  // maxValue/2
  bool binFilled = true;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (!densityMap.contains({rhmBin + 1, tMaxBin})) {
      binFilled = false;
      break;
    }
    rhmBin += 1;
    gridValue = densityMap.scaled({rhmBin, tMaxBin});
  }

  // Use linear approximation to find better z value for FWHM between bins
  double rightDensity = 0;
  if (binFilled) {
    rightDensity = densityMap.scaled({rhmBin, tMaxBin});
  }
  double leftDensity = densityMap.scaled({rhmBin - 1, tMaxBin});
  double deltaZ1 = m_cfg.spatialBinExtent * (maxValue / 2 - leftDensity) /
                   (rightDensity - leftDensity);

  int lhmBin = zMaxBin;
  gridValue = maxValue;
  binFilled = true;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (!densityMap.contains({lhmBin - 1, tMaxBin})) {
      binFilled = false;
      break;
    }
    lhmBin -= 1;
    gridValue = densityMap.scaled({lhmBin, tMaxBin});
  }

  // Use linear approximation to find better z value for FWHM between bins
  rightDensity = densityMap.scaled({lhmBin + 1, tMaxBin});
  if (binFilled) {
    leftDensity = densityMap.scaled({lhmBin, tMaxBin});
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
  densityMap.scaled(binFirstMax) = 0;
  auto secondMax = highestDensityEntry(densityMap);
  Bin binSecondMax = secondMax->first;
  double valueSecondMax = secondMax->second;
  double secondSum = 0;
  if (valueFirstMax - valueSecondMax < densityDeviation) {
    secondSum = getDensitySum(densityMap, binSecondMax);
  } else {
    // If the second maximum is not sufficiently large the third maximum won't
    // be either
    densityMap.scaled(binFirstMax) = valueFirstMax;
    return binFirstMax;
  }

  // Get the third highest maximum
  densityMap.scaled(binSecondMax) = 0;
  auto thirdMax = highestDensityEntry(densityMap);
  Bin binThirdMax = thirdMax->first;
  double valueThirdMax = thirdMax->second;
  double thirdSum = 0;
  if (valueFirstMax - valueThirdMax < densityDeviation) {
    thirdSum = getDensitySum(densityMap, binThirdMax);
  }

  // Revert back to original values
  densityMap.scaled(binFirstMax) = valueFirstMax;
  densityMap.scaled(binSecondMax) = valueSecondMax;

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
  // Add density from the bin.
  // Check if neighboring bins are part of the densityMap and add them (if they
  // are not part of the map, we assume them to be 0).
  return densityMap.scaled(bin) +
         densityMap.scaled({bin.first, bin.second - 1}) +
         densityMap.scaled({bin.first, bin.second + 1});
}

}  // namespace Acts
