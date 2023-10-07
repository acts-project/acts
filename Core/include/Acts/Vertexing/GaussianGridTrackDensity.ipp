// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>

template <int mainGridSize, int trkGridSize>
Acts::Result<float>
Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::getMaxZPosition(
    MainGridVector& mainGrid) const {
  if (mainGrid == MainGridVector::Zero()) {
    return VertexingError::EmptyInput;
  }

  int zbin = -1;
  if (!m_cfg.useHighestSumZPosition) {
    // Get bin with maximum content
    mainGrid.maxCoeff(&zbin);
  } else {
    // Get z position with highest density sum
    // of surrounding bins
    zbin = getHighestSumZPosition(mainGrid);
  }

  // Derive corresponding z value
  return (zbin - mainGridSize / 2 + 0.5f) * m_cfg.binSize;
}

template <int mainGridSize, int trkGridSize>
Acts::Result<std::pair<float, float>> Acts::GaussianGridTrackDensity<
    mainGridSize, trkGridSize>::getMaxZPositionAndWidth(MainGridVector&
                                                            mainGrid) const {
  // Get z maximum value
  auto maxZRes = getMaxZPosition(mainGrid);
  if (not maxZRes.ok()) {
    return maxZRes.error();
  }
  float maxZ = *maxZRes;

  // Get seed width estimate
  auto widthRes = estimateSeedWidth(mainGrid, maxZ);
  if (not widthRes.ok()) {
    return widthRes.error();
  }
  float width = *widthRes;
  std::pair<float, float> returnPair{maxZ, width};
  return returnPair;
}

template <int mainGridSize, int trkGridSize>
std::pair<int, typename Acts::GaussianGridTrackDensity<
                   mainGridSize, trkGridSize>::TrackGridVector>
Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::addTrack(
    const Acts::BoundTrackParameters& trk, MainGridVector& mainGrid) const {
  SymMatrix2 cov = trk.covariance().value().block<2, 2>(0, 0);
  float d0 = trk.parameters()[0];
  float z0 = trk.parameters()[1];

  // Calculate offset in d direction to central bin at z-axis
  int dOffset = static_cast<int>(std::floor(d0 / m_cfg.binSize - 0.5) + 1);
  // Calculate bin in z
  int zBin = int(z0 / m_cfg.binSize + mainGridSize / 2.);

  if (zBin < 0 || zBin >= mainGridSize) {
    return {-1, TrackGridVector::Zero()};
  }
  // Calculate the positions of the bin centers
  float binCtrD = dOffset * m_cfg.binSize;
  float binCtrZ = (zBin + 0.5f) * m_cfg.binSize - m_cfg.zMinMax;

  // Calculate the distance between IP values and their
  // corresponding bin centers
  float distCtrD = d0 - binCtrD;
  float distCtrZ = z0 - binCtrZ;

  // Check if current track does affect grid density
  // in central bins at z-axis
  if (std::abs(dOffset) > (trkGridSize - 1) / 2.) {
    // Current track is too far away to contribute
    // to track density at z-axis bins
    return {-1, TrackGridVector::Zero()};
  }

  // Create the track grid
  TrackGridVector trackGrid = createTrackGrid(dOffset, cov, distCtrD, distCtrZ);
  // Add the track grid to the main grid
  addTrackGridToMainGrid(zBin, trackGrid, mainGrid);

  return {zBin, trackGrid};
}

template <int mainGridSize, int trkGridSize>
void Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::
    addTrackGridToMainGrid(int zBin, const TrackGridVector& trkGrid,
                           MainGridVector& mainGrid) const {
  modifyMainGridWithTrackGrid(zBin, trkGrid, mainGrid, +1);
}

template <int mainGridSize, int trkGridSize>
void Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::
    removeTrackGridFromMainGrid(int zBin, const TrackGridVector& trkGrid,
                                MainGridVector& mainGrid) const {
  modifyMainGridWithTrackGrid(zBin, trkGrid, mainGrid, -1);
}

template <int mainGridSize, int trkGridSize>
void Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::
    modifyMainGridWithTrackGrid(int zBin, const TrackGridVector& trkGrid,
                                MainGridVector& mainGrid,
                                int modifyModeSign) const {
  int width = (trkGridSize - 1) / 2;
  // Overlap left
  int leftOL = zBin - width;
  // Overlap right
  int rightOL = zBin + width - mainGridSize + 1;
  if (leftOL < 0) {
    int totalTrkSize = trkGridSize + leftOL;
    mainGrid.segment(0, totalTrkSize) +=
        modifyModeSign * trkGrid.segment(-leftOL, totalTrkSize);
    return;
  }
  if (rightOL > 0) {
    int totalTrkSize = trkGridSize - rightOL;
    mainGrid.segment(mainGridSize - totalTrkSize, totalTrkSize) +=
        modifyModeSign * trkGrid.segment(0, totalTrkSize);
    return;
  }

  mainGrid.segment(zBin - width, trkGridSize) += modifyModeSign * trkGrid;
}

template <int mainGridSize, int trkGridSize>
typename Acts::GaussianGridTrackDensity<mainGridSize,
                                        trkGridSize>::TrackGridVector
Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::createTrackGrid(
    int offset, const Acts::SymMatrix2& cov, float distCtrD,
    float distCtrZ) const {
  TrackGridVector trackGrid(TrackGridVector::Zero());

  int i = (trkGridSize - 1) / 2 + offset;
  float d = (i - static_cast<float>(trkGridSize) / 2 + 0.5f) * m_cfg.binSize;

  // Loop over columns
  for (int j = 0; j < trkGridSize; j++) {
    float z = (j - static_cast<float>(trkGridSize) / 2 + 0.5f) * m_cfg.binSize;
    trackGrid(j) = normal2D(d + distCtrD, z + distCtrZ, cov);
  }
  return trackGrid;
}

template <int mainGridSize, int trkGridSize>
Acts::Result<float>
Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::estimateSeedWidth(
    MainGridVector& mainGrid, float maxZ) const {
  if (mainGrid == MainGridVector::Zero()) {
    return VertexingError::EmptyInput;
  }
  // Get z bin of max density z value
  int zBin = int(maxZ / m_cfg.binSize + mainGridSize / 2.);

  const float maxValue = mainGrid(zBin);
  float gridValue = mainGrid(zBin);

  // Find right half-maximum bin
  int rhmBin = zBin;
  while (gridValue > maxValue / 2) {
    rhmBin += 1;
    gridValue = mainGrid(rhmBin);
  }

  // Use linear approximation to find better z value for FWHM between bins
  float deltaZ1 = (maxValue / 2 - mainGrid(rhmBin - 1)) *
                  (m_cfg.binSize / (mainGrid(rhmBin - 1) - mainGrid(rhmBin)));

  // Find left half-maximum bin
  int lhmBin = zBin;
  gridValue = mainGrid(zBin);
  while (gridValue > maxValue / 2) {
    lhmBin -= 1;
    gridValue = mainGrid(lhmBin);
  }

  // Use linear approximation to find better z value for FWHM between bins
  float deltaZ2 = (maxValue / 2 - mainGrid(lhmBin + 1)) *
                  (m_cfg.binSize / (mainGrid(rhmBin + 1) - mainGrid(rhmBin)));

  // Approximate FWHM
  float fwhm =
      rhmBin * m_cfg.binSize - deltaZ1 - lhmBin * m_cfg.binSize - deltaZ2;

  // FWHM = 2.355 * sigma
  float width = fwhm / 2.355f;

  return std::isnormal(width) ? width : 0.0f;
}

template <int mainGridSize, int trkGridSize>
float Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::normal2D(
    float d, float z, const Acts::SymMatrix2& cov) const {
  float det = cov.determinant();
  float coef = 1 / (2 * M_PI * std::sqrt(det));
  float expo =
      -1 / (2 * det) *
      (cov(1, 1) * d * d - d * z * (cov(0, 1) + cov(1, 0)) + cov(0, 0) * z * z);
  return coef * std::exp(expo);
}

template <int mainGridSize, int trkGridSize>
int Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::
    getHighestSumZPosition(MainGridVector& mainGrid) const {
  // Checks the first (up to) 3 density maxima, if they are close, checks which
  // one has the highest surrounding density sum (the two neighboring bins)

  // The global maximum
  int zbin = -1;
  mainGrid.maxCoeff(&zbin);
  int zFirstMax = zbin;
  double firstDensity = mainGrid(zFirstMax);
  double firstSum = getDensitySum(mainGrid, zFirstMax);

  // Get the second highest maximum
  mainGrid[zFirstMax] = 0;
  mainGrid.maxCoeff(&zbin);
  int zSecondMax = zbin;
  double secondDensity = mainGrid(zSecondMax);
  double secondSum = 0;
  if (firstDensity - secondDensity <
      firstDensity * m_cfg.maxRelativeDensityDev) {
    secondSum = getDensitySum(mainGrid, zSecondMax);
  }

  // Get the third highest maximum
  mainGrid[zSecondMax] = 0;
  mainGrid.maxCoeff(&zbin);
  int zThirdMax = zbin;
  double thirdDensity = mainGrid(zThirdMax);
  double thirdSum = 0;
  if (firstDensity - thirdDensity <
      firstDensity * m_cfg.maxRelativeDensityDev) {
    thirdSum = getDensitySum(mainGrid, zThirdMax);
  }

  // Revert back to original values
  mainGrid[zFirstMax] = firstDensity;
  mainGrid[zSecondMax] = secondDensity;

  // Return the z-bin position of the highest density sum
  if (secondSum > firstSum && secondSum > thirdSum) {
    return zSecondMax;
  }
  if (thirdSum > secondSum && thirdSum > firstSum) {
    return zThirdMax;
  }
  return zFirstMax;
}

template <int mainGridSize, int trkGridSize>
double Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::getDensitySum(
    const MainGridVector& mainGrid, int pos) const {
  double sum = mainGrid(pos);
  // Sum up only the density contributions from the
  // neighboring bins if they are still within bounds
  if (pos - 1 >= 0) {
    sum += mainGrid(pos - 1);
  }
  if (pos + 1 <= mainGrid.size() - 1) {
    sum += mainGrid(pos + 1);
  }
  return sum;
}
