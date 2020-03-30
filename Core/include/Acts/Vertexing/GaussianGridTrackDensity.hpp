// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <iostream>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

template <int nDbins = 80, int nZbins = 2000, int nTrkDbins = 7,
          int nTrkZbins = 7>
class GaussianGridTrackDensity {
 public:
  struct Config {
    Config(double zMinMax = 100, double dMinMax = 4)
        : zMinMax(zMinMax), dMinMax(dMinMax) {
      dBinSize = 2. * dMinMax / nDbins;
      zBinSize = 2. * zMinMax / nZbins;
    }

    // Min and max z value of big grid
    double zMinMax;  // mm
    // Min and max d value of big grid
    double dMinMax;  // mm

    // Z size of one single bin in grid
    double zBinSize;  // mm
    // Z size of one single bin in grid
    double dBinSize;  // mm
  };

  GaussianGridTrackDensity(const Config& cfg) : m_cfg(cfg) {}

  void addTrackToGrid(const Acts::BoundParameters& trk) const;

 private:
  ActsMatrixF<nTrkDbins, nTrkZbins> createTrackGrid(
      const ActsSymMatrixD<2>& cov, double distCtrD, double distCtrZ) const;

  double normal2D(double d, double z, const ActsSymMatrixD<2>& cov) const;

  Config m_cfg;
};

}  // namespace Acts

#include "Acts/Vertexing/GaussianGridTrackDensity.ipp"