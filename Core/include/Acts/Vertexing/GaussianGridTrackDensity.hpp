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

template <int mainGridSize = 2000, int trkGridSize = 15>
class GaussianGridTrackDensity {
  // Assert odd trkGridSize
  static_assert(trkGridSize % 2);
  // Assert bigger main grid than track grid
  static_assert(mainGridSize > trkGridSize);

 public:
  struct Config {
    Config(double zMinMax = 100) : zMinMax(zMinMax) {
      binSize = 2. * zMinMax / mainGridSize;
    }

    // Min and max z value of big grid
    double zMinMax;  // mm

    // Z size of one single bin in grid
    double binSize;  // mm
  };

  GaussianGridTrackDensity(const Config& cfg) : m_cfg(cfg) {}

  void addTrack(const Acts::BoundParameters& trk,
                Acts::ActsVectorF<mainGridSize>& mainGrid) const;

 private:
  void addTrackGridToMainGrid(int zBin, const ActsVectorF<trkGridSize>& trkGrid,
                              Acts::ActsVectorF<mainGridSize>& mainGrid) const;

  ActsVectorF<trkGridSize> createTrackGrid(int offset,
                                           const ActsSymMatrixD<2>& cov,
                                           double distCtrD,
                                           double distCtrZ) const;

  double normal2D(double d, double z, const ActsSymMatrixD<2>& cov) const;

  Config m_cfg;
};

}  // namespace Acts

#include "Acts/Vertexing/GaussianGridTrackDensity.ipp"