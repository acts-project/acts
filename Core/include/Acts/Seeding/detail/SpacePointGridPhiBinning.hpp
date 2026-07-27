// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

namespace Acts::detail {

/// Parameters to derive the number of phi bins of a space point grid from the
/// minimum track transverse momentum and the magnetic field.
struct SpacePointGridPhiBinningConfig {
  /// minimum pT
  float minPt = 0;
  /// magnetic field
  float bFieldInZ = 0;
  /// maximum extension of sensitive detector layer relevant for seeding as
  /// distance from x=y=0 (i.e. in r)
  float rMax = 0;
  /// maximum distance in r from middle space point to bottom or top
  /// space point
  float deltaRMax = 0;
  /// maximum impact parameter
  float impactMax = 0;
  /// Multiplicator for the number of phi-bins. The minimum number of phi-bins
  /// depends on min_pt, magnetic field: 2*pi/(minPT particle phi-deflection).
  /// phiBinDeflectionCoverage is a multiplier for this number. If
  /// numPhiNeighbors (in the configuration of the BinFinders) is configured
  /// to return 1 neighbor on either side of the current phi-bin (and you want
  /// to cover the full phi-range of minPT), leave this at 1.
  int phiBinDeflectionCoverage = 1;
  /// maximum number of phi bins
  int maxPhiBins = 10000;
};

/// Compute the number of phi bins of a space point grid such that each bin
/// covers the maximum expected azimuthal deflection of a minimum-pT track
/// between the innermost and outermost radius relevant for seeding, including
/// the effect of the maximum impact parameter.
/// @param config Parameters the phi bin count is derived from
/// @param logger Logger instance for debugging output
/// @return The number of phi bins, capped at `config.maxPhiBins`
int computeSpacePointGridPhiBins(const SpacePointGridPhiBinningConfig& config,
                                 const Logger& logger = getDummyLogger());

}  // namespace Acts::detail
