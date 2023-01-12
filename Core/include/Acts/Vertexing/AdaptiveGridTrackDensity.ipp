// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>

template <int trkGridSize>
Acts::Result<float>
Acts::AdaptiveGridTrackDensity<trkGridSize>::getMaxZPosition(
    std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const {
  if (mainGridDensity.empty() || mainGridZValues.empty()) {
    return VertexingError::EmptyInput;
  }

  int zGridPos = -1;
  if (!m_cfg.useHighestSumZPosition) {
    zGridPos = std::distance(
        mainGridDensity.begin(),
        std::max_element(mainGridDensity.begin(), mainGridDensity.end()));
  } else {
    // Get z position with highest density sum
    // of surrounding bins
    zGridPos = getHighestSumZPosition(mainGridDensity, mainGridZValues);
  }

  int zbin = mainGridZValues[zGridPos];
  // Derive corresponding z value
  int sign = (zbin > 0) ? +1 : -1;
  return (zbin + sign * 0.5f) * m_cfg.binSize;
}

template <int trkGridSize>
Acts::Result<std::pair<float, float>>
Acts::AdaptiveGridTrackDensity<trkGridSize>::getMaxZPositionAndWidth(
    std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const {
  // Get z maximum value
  auto maxZRes = getMaxZPosition(mainGridDensity, mainGridZValues);
  if (not maxZRes.ok()) {
    return maxZRes.error();
  }
  float maxZ = *maxZRes;

  // Get seed width estimate
  auto widthRes = estimateSeedWidth(mainGridDensity, mainGridZValues, maxZ);
  if (not widthRes.ok()) {
    return widthRes.error();
  }
  float width = *widthRes;
  std::pair<float, float> returnPair{maxZ, width};
  return returnPair;
}

template <int trkGridSize>
std::pair<int,
          typename Acts::AdaptiveGridTrackDensity<trkGridSize>::TrackGridVector>
Acts::AdaptiveGridTrackDensity<trkGridSize>::addTrack(
    const Acts::BoundTrackParameters& trk, std::vector<float>& mainGridDensity,
    std::vector<int>& mainGridZValues) const {
  SymMatrix2 cov = trk.covariance().value().block<2, 2>(0, 0);
  float d0 = trk.parameters()[0];
  float z0 = trk.parameters()[1];

  // Calculate offset in d direction to central bin at z-axis
  int dOffset = std::floor(d0 / m_cfg.binSize - 0.5) + 1;
  // Calculate bin in z
  int zBin = int(z0 / m_cfg.binSize);

  // Calculate the positions of the bin centers
  float binCtrD = dOffset * m_cfg.binSize;

  int sign = (z0 > 0) ? +1 : -1;
  float binCtrZ = (zBin + sign * 0.5f) * m_cfg.binSize;

  // Calculate the distance between IP values and their
  // corresponding bin centers
  float distCtrD = d0 - binCtrD;
  float distCtrZ = z0 - binCtrZ;

  TrackGridVector trackGrid(TrackGridVector::Zero());

  // Check if current track does affect grid density
  // in central bins at z-axis
  if (std::abs(dOffset) > (trkGridSize - 1) / 2.) {
    return {0, trackGrid};
  }

  // Create the track grid
  trackGrid = createTrackGrid(dOffset, cov, distCtrD, distCtrZ);

  std::vector<int> zBinValues;

  int startEnd = int(trkGridSize - 1) / 2;

  for (int i = 0; i < trkGridSize; i++) {
    zBinValues.push_back(int(zBin + (i - startEnd)));
  }

  for (int i = 0; i < trkGridSize; i++) {
    int z = zBinValues[i];

    // Check if track density already exists at current z position
    auto findIter =
        std::find(mainGridZValues.begin(), mainGridZValues.end(), z);

    if (findIter != mainGridZValues.end()) {
      // Z bin already exists
      mainGridDensity[std::distance(mainGridZValues.begin(), findIter)] +=
          trackGrid[i];
    } else {
      // Create new z bin
      auto it =
          std::upper_bound(mainGridZValues.begin(), mainGridZValues.end(), z);
      mainGridDensity.insert(
          mainGridDensity.begin() + std::distance(mainGridZValues.begin(), it),
          trackGrid[i]);
      mainGridZValues.insert(it, z);
    }
  }

  return {zBin, trackGrid};
}

template <int trkGridSize>
void Acts::AdaptiveGridTrackDensity<trkGridSize>::removeTrackGridFromMainGrid(
    int zBin, const TrackGridVector& trkGrid,
    std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const {
  // Find position of current z bin in mainGridZValues
  auto findIter =
      std::find(mainGridZValues.begin(), mainGridZValues.end(), zBin);
  // Calculate corresponding index in mainGridDensity
  int densityIdx = std::distance(mainGridZValues.begin(), findIter);

  // Go over trkGrid and remove it from mainDensityGrid
  int startEnd = int((trkGridSize - 1) / 2);
  for (int i = 0; i < trkGridSize; i++) {
    mainGridDensity[int(densityIdx + (i - startEnd))] -= trkGrid[i];
  }
}

template <int trkGridSize>
typename Acts::AdaptiveGridTrackDensity<trkGridSize>::TrackGridVector
Acts::AdaptiveGridTrackDensity<trkGridSize>::createTrackGrid(
    int offset, const Acts::SymMatrix2& cov, float distCtrD,
    float distCtrZ) const {
  TrackGridVector trackGrid(TrackGridVector::Zero());

  float i = (trkGridSize - 1) / 2 + offset;
  float d = (i - static_cast<float>(trkGridSize) / 2 + 0.5f) * m_cfg.binSize;

  // Loop over columns
  for (int j = 0; j < trkGridSize; j++) {
    float z = (j - static_cast<float>(trkGridSize) / 2 + 0.5f) * m_cfg.binSize;
    trackGrid(j) = normal2D(d + distCtrD, z + distCtrZ, cov);
  }
  return trackGrid;
}

template <int trkGridSize>
Acts::Result<float>
Acts::AdaptiveGridTrackDensity<trkGridSize>::estimateSeedWidth(
    const std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues, float maxZ) const {
  if (mainGridDensity.empty() || mainGridZValues.empty()) {
    return VertexingError::EmptyInput;
  }
  // Get z bin of max density z value
  int sign = (maxZ > 0) ? +1 : -1;
  int zMaxGridBin = int(maxZ / m_cfg.binSize - sign * 0.5f);

  // Find location in mainGridZValues
  auto findIter =
      std::find(mainGridZValues.begin(), mainGridZValues.end(), zMaxGridBin);
  int zBin = std::distance(mainGridZValues.begin(), findIter);

  const float maxValue = mainGridDensity[zBin];
  float gridValue = mainGridDensity[zBin];

  // Find right half-maximum bin
  int rhmBin = zBin;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continous z values
    if ((zMaxGridBin + (rhmBin - zBin)) != mainGridZValues[rhmBin]) {
      break;
    }
    rhmBin += 1;
    if (rhmBin == int(mainGridDensity.size())) {
      break;
    }
    gridValue = mainGridDensity[rhmBin];
  }

  // Use linear approximation to find better z value for FWHM between bins
  float deltaZ1 =
      (maxValue / 2 - mainGridDensity[rhmBin - 1]) *
      (m_cfg.binSize / (mainGridDensity[rhmBin - 1] - mainGridDensity[rhmBin]));
  // Find left half-maximum bin
  int lhmBin = zBin;
  gridValue = mainGridDensity[zBin];
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continous z values
    if ((zMaxGridBin + (lhmBin - zBin)) != mainGridZValues[lhmBin]) {
      break;
    }
    lhmBin -= 1;
    if (lhmBin < 0) {
      break;
    }
    gridValue = mainGridDensity[lhmBin];
  }

  // Use linear approximation to find better z value for FWHM between bins
  float deltaZ2 =
      (maxValue / 2 - mainGridDensity[lhmBin + 1]) *
      (m_cfg.binSize / (mainGridDensity[rhmBin + 1] - mainGridDensity[rhmBin]));

  // Approximate FWHM
  float fwhm =
      rhmBin * m_cfg.binSize - deltaZ1 - lhmBin * m_cfg.binSize - deltaZ2;

  // FWHM = 2.355 * sigma
  float width = fwhm / 2.355f;

  return std::isnormal(width) ? width : 0.0f;
}

template <int trkGridSize>
float Acts::AdaptiveGridTrackDensity<trkGridSize>::normal2D(
    float d, float z, const Acts::SymMatrix2& cov) const {
  float det = cov.determinant();
  float coef = 1 / (2 * M_PI * std::sqrt(det));
  float expo =
      -1 / (2 * det) *
      (cov(1, 1) * d * d - d * z * (cov(0, 1) + cov(1, 0)) + cov(0, 0) * z * z);
  return coef * std::exp(expo);
}

template <int trkGridSize>
int Acts::AdaptiveGridTrackDensity<trkGridSize>::getHighestSumZPosition(
    std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const {
  // Checks the first (up to) 3 density maxima, if they are close, checks which
  // one has the highest surrounding density sum (the two neighboring bins)

  // The global maximum
  unsigned int zGridPos = std::distance(
      mainGridDensity.begin(),
      std::max_element(mainGridDensity.begin(), mainGridDensity.end()));

  unsigned int zFirstMax = zGridPos;
  double firstDensity = mainGridDensity[zFirstMax];
  double firstSum = getDensitySum(mainGridDensity, mainGridZValues, zFirstMax);

  // Get the second highest maximum
  mainGridDensity[zFirstMax] = 0;
  zGridPos = std::distance(
      mainGridDensity.begin(),
      std::max_element(mainGridDensity.begin(), mainGridDensity.end()));
  unsigned int zSecondMax = zGridPos;
  double secondDensity = mainGridDensity[zSecondMax];
  double secondSum = 0;
  if (firstDensity - secondDensity <
      firstDensity * m_cfg.maxRelativeDensityDev) {
    secondSum = getDensitySum(mainGridDensity, mainGridZValues, zSecondMax);
  }

  // Get the third highest maximum
  mainGridDensity[zSecondMax] = 0;
  zGridPos = std::distance(
      mainGridDensity.begin(),
      std::max_element(mainGridDensity.begin(), mainGridDensity.end()));
  unsigned int zThirdMax = zGridPos;
  double thirdDensity = mainGridDensity[zThirdMax];
  double thirdSum = 0;
  if (firstDensity - thirdDensity <
      firstDensity * m_cfg.maxRelativeDensityDev) {
    thirdSum = getDensitySum(mainGridDensity, mainGridZValues, zThirdMax);
  }

  // Revert back to original values
  mainGridDensity[zFirstMax] = firstDensity;
  mainGridDensity[zSecondMax] = secondDensity;

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
double Acts::AdaptiveGridTrackDensity<trkGridSize>::getDensitySum(
    const std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues, unsigned int pos) const {
  double sum = mainGridDensity[pos];
  // Sum up only the density contributions from the
  // neighboring bins if they are still within bounds
  if (0 < pos) {
    // Check if we are still operating on continous z values
    if (mainGridZValues[pos] - mainGridZValues[pos - 1] == 1) {
      sum += mainGridDensity[pos - 1];
    }
  }
  if (pos + 1 <= mainGridDensity.size() - 1) {
    // Check if we are still operating on continous z values
    if (mainGridZValues[pos + 1] - mainGridZValues[pos] == 1) {
      sum += mainGridDensity[pos + 1];
    }
  }
  return sum;
}
