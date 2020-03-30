// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <int nDbins, int nZbins, int nTrkDbins, int nTrkZbins>
void Acts::GaussianGridTrackDensity<nDbins, nZbins, nTrkDbins, nTrkZbins>::
    addTrackToGrid(const Acts::BoundParameters& trk) const {
  ActsSymMatrixD<2> cov = trk.covariance()->block<2, 2>(0, 0);
  double d0 = trk.parameters()[0];
  double z0 = trk.parameters()[1];

  std::cout << "IP: " << d0 << ", " << z0 << std::endl;

  // Get the bin numbers corresponding to d0/z0 values
  int d0Bin = int(d0 / m_cfg.dBinSize + nDbins / 2);
  int z0Bin = int(z0 / m_cfg.zBinSize + nZbins / 2);

  std::cout << "bins: " << d0Bin << ", " << z0Bin << std::endl;

  std::cout << m_cfg.dBinSize << ", " << m_cfg.zBinSize << std::endl;
  std::cout << m_cfg.dMinMax << ", " << m_cfg.zMinMax << std::endl;

  // Find the bin center values
  double binCtrD = (d0Bin + 0.5) * m_cfg.dBinSize - m_cfg.dMinMax;
  double binCtrZ = (z0Bin + 0.5) * m_cfg.zBinSize - m_cfg.zMinMax;

  std::cout << "binCtr: " << binCtrD << ", " << binCtrZ << std::endl;

  // Compute the distance between d0/z0 and the bin centers
  double distCtrD = d0 - binCtrD;
  double distCtrZ = z0 - binCtrZ;

  // Create the track grid
  ActsMatrixF<nTrkDbins, nTrkZbins> trackGrid =
      createTrackGrid(cov, distCtrD, distCtrZ);

  std::cout << "dist: " << distCtrD << ", " << distCtrZ << std::endl;
  std::cout << trackGrid << std::endl;
}

template <int nDbins, int nZbins, int nTrkDbins, int nTrkZbins>
Acts::ActsMatrixF<nTrkDbins, nTrkZbins> Acts::GaussianGridTrackDensity<
    nDbins, nZbins, nTrkDbins,
    nTrkZbins>::createTrackGrid(const Acts::ActsSymMatrixD<2>& cov,
                                double distCtrD, double distCtrZ) const {
  ActsMatrixF<nTrkDbins, nTrkZbins> trackGrid(
      ActsMatrixF<nTrkDbins, nTrkZbins>::Zero());
  // Loop over columns
  for (int j = 0; j < nTrkZbins; j++) {
    // Loop over rows
    for (int i = 0; i < nTrkDbins; i++) {
      double d = (i - (double)nTrkDbins / 2. + 0.5) * m_cfg.dBinSize;
      double z = (j - (double)nTrkZbins / 2. + 0.5) * m_cfg.zBinSize;
      trackGrid(i, j) = normal2D(d + distCtrD, z + distCtrZ, cov);
    }
  }
  return trackGrid;
}

template <int nDbins, int nZbins, int nTrkDbins, int nTrkZbins>
double
Acts::GaussianGridTrackDensity<nDbins, nZbins, nTrkDbins, nTrkZbins>::normal2D(
    double d, double z, const Acts::ActsSymMatrixD<2>& cov) const {
  double det = cov.determinant();
  double coef = 1. / (2. * M_PI * std::sqrt(det));
  double expo =
      -1. / (2. * det) *
      (cov(1, 1) * d * d - d * z * (cov(0, 1) + cov(1, 0)) + cov(0, 0) * z * z);
  return coef * std::exp(expo);
}