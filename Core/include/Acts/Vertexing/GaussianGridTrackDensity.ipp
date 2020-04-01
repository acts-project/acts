// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "Acts/Vertexing/VertexingError.hpp"

template <int mainGridSize, int trkGridSize>
Acts::Result<float>
Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::getMaxZPosition(
    const Acts::ActsVectorF<mainGridSize>& mainGrid) const {
  if (mainGrid == ActsVectorF<mainGridSize>::Zero()) {
    return VertexingError::EmptyInput;
  }

  // Get bin with maximum content
  int zbin = -1;
  mainGrid.maxCoeff(&zbin);
  // Derive corresponding z value
  return (zbin - mainGridSize / 2) * m_cfg.binSize;
}

template <int mainGridSize, int trkGridSize>
std::pair<int, Acts::ActsVectorF<trkGridSize>>
Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::addTrack(
    const Acts::BoundParameters& trk,
    Acts::ActsVectorF<mainGridSize>& mainGrid) const {
  ActsSymMatrixD<2> cov = trk.covariance()->block<2, 2>(0, 0);
  float d0 = trk.parameters()[0];
  float z0 = trk.parameters()[1];

  // Calculate offset in d direction to central bin at z = 0
  int dOffset = std::floor(d0 / m_cfg.binSize - 0.5) + 1;
  // Calculate bin in z
  int zBin = int(z0 / m_cfg.binSize + mainGridSize / 2.);

  if (zBin < 0 || zBin >= mainGridSize) {
    return {-1, ActsVectorF<trkGridSize>::Zero()};
  }
  // Calculate the positions of the bin centers
  float binCtrD = dOffset * m_cfg.binSize;
  float binCtrZ = (zBin + 0.5) * m_cfg.binSize - m_cfg.zMinMax;

  // Calculate the distance between IP values and their
  // corresponding bin centers
  float distCtrD = d0 - binCtrD;
  float distCtrZ = z0 - binCtrZ;

  // Check if current track does affect grid density
  // in central bins at z = 0
  if ((std::abs(dOffset) > trkGridSize - 1) / 2.) {
    // Current track is too far away to contribute
    // to track density at z = 0 bins
    return {-1, ActsVectorF<trkGridSize>::Zero()};
  }

  // Create the track grid
  ActsVectorF<trkGridSize> trackGrid =
      createTrackGrid(dOffset, cov, distCtrD, distCtrZ);
  // Add the track grid to the main grid
  addTrackGridToMainGrid(zBin, trackGrid, mainGrid);

  return {zBin, trackGrid};
}

template <int mainGridSize, int trkGridSize>
void Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::
    addTrackGridToMainGrid(int zBin,
                           const Acts::ActsVectorF<trkGridSize>& trkGrid,
                           Acts::ActsVectorF<mainGridSize>& mainGrid) const {
  modifyMainGridWithTrackGrid(zBin, trkGrid, mainGrid, +1);
}

template <int mainGridSize, int trkGridSize>
void Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::
    removeTrackGridFromMainGrid(
        int zBin, const Acts::ActsVectorF<trkGridSize>& trkGrid,
        Acts::ActsVectorF<mainGridSize>& mainGrid) const {
  modifyMainGridWithTrackGrid(zBin, trkGrid, mainGrid, -1);
}

template <int mainGridSize, int trkGridSize>
void Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::
    modifyMainGridWithTrackGrid(int zBin,
                                const Acts::ActsVectorF<trkGridSize>& trkGrid,
                                Acts::ActsVectorF<mainGridSize>& mainGrid,
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
Acts::ActsVectorF<trkGridSize>
Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::createTrackGrid(
    int offset, const Acts::ActsSymMatrixD<2>& cov, float distCtrD,
    float distCtrZ) const {
  ActsVectorF<trkGridSize> trackGrid(ActsVectorF<trkGridSize>::Zero());

  int i = (trkGridSize - 1) / 2. + offset;
  float d = (i - (float)trkGridSize / 2. + 0.5) * m_cfg.binSize;

  // Loop over columns
  for (int j = 0; j < trkGridSize; j++) {
    float z = (j - (float)trkGridSize / 2. + 0.5) * m_cfg.binSize;
    trackGrid(j) = normal2D(d + distCtrD, z + distCtrZ, cov);
  }
  return trackGrid;
}

template <int mainGridSize, int trkGridSize>
float Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::normal2D(
    float d, float z, const Acts::ActsSymMatrixD<2>& cov) const {
  float det = cov.determinant();
  float coef = 1. / (2. * M_PI * std::sqrt(det));
  float expo =
      -1. / (2. * det) *
      (cov(1, 1) * d * d - d * z * (cov(0, 1) + cov(1, 0)) + cov(0, 0) * z * z);
  return coef * std::exp(expo);
}