// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <int mainGridSize, int trkGridSize>
void Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::addTrack(
    const Acts::BoundParameters& trk,
    Acts::ActsVectorF<mainGridSize>& mainGrid) const {
  ActsSymMatrixD<2> cov = trk.covariance()->block<2, 2>(0, 0);
  double d0 = trk.parameters()[0];
  double z0 = trk.parameters()[1];

  std::cout << "IP: " << d0 << ", " << z0 << std::endl;

  // Calculate offset in d direction to central bin at z = 0
  int dOffset = std::floor(d0 / m_cfg.binSize - 0.5) + 1;

  std::cout << "offset: " << dOffset << std::endl;

  // Calculate bin in z
  int zBin = int(z0 / m_cfg.binSize + mainGridSize / 2.);

  // Calculate the positions of the bin centers
  double binCtrD = dOffset * m_cfg.binSize;
  double binCtrZ = (zBin + 0.5) * m_cfg.binSize - m_cfg.zMinMax;

  // Calculate the distance between IP values and their
  // corresponding bin centers
  double distCtrD = d0 - binCtrD;
  double distCtrZ = z0 - binCtrZ;

  // Check if current track does affect grid density
  // in central bins at z = 0
  if ((std::abs(dOffset) > trkGridSize - 1) / 2.) {
    // Current track is too far away to contribute
    // to track density at z = 0 bins
    return;
  }

  // Create the track grid
  ActsVectorF<trkGridSize> trackGrid =
      createTrackGrid(dOffset, cov, distCtrD, distCtrZ);
  // TODO: test correct output.. symmetry?
  std::cout << trackGrid << std::endl;

  std::cout << "mainGrid before: " << mainGrid << std::endl;
  addTrackGridToMainGrid(zBin, trackGrid, mainGrid);
  std::cout << "mainGrid after: " << mainGrid << std::endl;
}

template <int mainGridSize, int trkGridSize>
void Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::
    addTrackGridToMainGrid(int zBin,
                           const Acts::ActsVectorF<trkGridSize>& trkGrid,
                           Acts::ActsVectorF<mainGridSize>& mainGrid) const {
  int width = (trkGridSize - 1) / 2;
  // Overlap left
  int leftOL = zBin - width;
  // Overlap right
  int rightOL = zBin + width - mainGridSize;
  if (leftOL < 0) {
    int totalTrkSize = trkGridSize + leftOL;
    mainGrid.segment(0, totalTrkSize) += trkGrid.segment(-leftOL, totalTrkSize);
    return;
  }
  if (rightOL > 0) {
    int totalTrkSize = trkGridSize - rightOL;
    mainGrid.segment(mainGridSize - totalTrkSize, totalTrkSize) +=
        trkGrid.segment(0, totalTrkSize);
    return;
  }

  mainGrid.segment(zBin - width, trkGridSize) += trkGrid;
}

template <int mainGridSize, int trkGridSize>
Acts::ActsVectorF<trkGridSize>
Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::createTrackGrid(
    int offset, const Acts::ActsSymMatrixD<2>& cov, double distCtrD,
    double distCtrZ) const {
  ActsVectorF<trkGridSize> trackGrid(ActsVectorF<trkGridSize>::Zero());

  int i = (trkGridSize - 1) / 2. - offset;
  double d = (i - (double)trkGridSize / 2. + 0.5) * m_cfg.binSize;

  // Loop over columns
  for (int j = 0; j < trkGridSize; j++) {
    double z = (j - (double)trkGridSize / 2. + 0.5) * m_cfg.binSize;
    trackGrid(j) = normal2D(d + distCtrD, z + distCtrZ, cov);
  }
  return trackGrid;
}

template <int mainGridSize, int trkGridSize>
double Acts::GaussianGridTrackDensity<mainGridSize, trkGridSize>::normal2D(
    double d, double z, const Acts::ActsSymMatrixD<2>& cov) const {
  double det = cov.determinant();
  double coef = 1. / (2. * M_PI * std::sqrt(det));
  double expo =
      -1. / (2. * det) *
      (cov(1, 1) * d * d - d * z * (cov(0, 1) + cov(1, 0)) + cov(0, 0) * z * z);
  return coef * std::exp(expo);
}