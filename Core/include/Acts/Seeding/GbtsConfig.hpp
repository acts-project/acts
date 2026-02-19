// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

#include <cstdint>
#include <string>

namespace Acts::Experimental {

/// Configuration options for the GBTs seed finder.
struct GbtsConfig {
  // GbtsSeedingAlgorithm options
  /// Enable beam spot correction.
  bool beamSpotCorrection = false;

  // Path to the connector configuration file that defines the layer connections
  /// Connector configuration file path.
  std::string connectorInputFile;

  /// Look-up table input file path.
  std::string lutInputFile;

  // SeedFinderGbts option
  /// Enable Large Radius Tracking mode.
  bool lrtMode = false;
  /// Use machine learning features (e.g., cluster width).
  bool useMl = false;
  /// Match seeds before creating them.
  bool matchBeforeCreate = false;
  /// Use legacy tuning parameters.
  bool useOldTunings = false;
  /// Tau ratio cut threshold.
  float tauRatioCut = 0.007;
  /// Tau ratio precut threshold.
  float tauRatioPrecut = 0.009f;
  /// Eta bin width override (0 uses default from connection file).
  float etaBinOverride =
      0.0f;  // specify non-zero to override eta bin width from connection file
             // (default 0.2 in createLinkingScheme.py)

  /// Maximum number of phi slices.
  float nMaxPhiSlice = 53;  // used to calculate phi slices
  /// Minimum transverse momentum.
  float minPt = 1.0f * UnitConstants::GeV;
  /// Phi slice width (derived in CreateSeeds function).
  float phiSliceWidth{};  // derived in CreatSeeds function

  // graph building options
  /// Transverse momentum coefficient (~0.3*B/2 - assumes nominal field of 2*T).
  double ptCoeff = 0.29997 * 1.9972 / 2.0;
  /// Use eta binning from geometry structure.
  bool useEtaBinning = true;
  /// Apply RZ cuts on doublets.
  bool doubletFilterRZ = true;
  /// Maximum number of Gbts edges/doublets.
  std::uint32_t nMaxEdges = 2000000;
  /// Minimum delta radius between layers.
  float minDeltaRadius = 2.0;

  // GbtsTrackingFilter options
  /// Multiple scattering sigma (for 900 MeV track at eta=0).
  float sigmaMS = 0.016;
  /// Radiation length fraction per layer (2.5% per layer).
  float radLen = 0.025;

  /// Measurement uncertainty in x direction.
  float sigmaX = 0.08;
  /// Measurement uncertainty in y direction.
  float sigmaY = 0.25;

  /// Measurement weight in x direction.
  float weightX = 0.5;
  /// Measurement weight in y direction.
  float weightY = 0.5;

  /// Maximum delta chi2 in x direction.
  float maxDChi2X = 5.0;
  /// Maximum delta chi2 in y direction.
  float maxDChi2Y = 6.0;

  /// Chi2 penalty for adding a hit.
  float addHit = 14.0;

  /// Maximum track curvature.
  float maxCurvature = 1e-3f;
  /// Maximum longitudinal impact parameter.
  float maxZ0 = 170.0;

  // Seed extraction options
  /// Minimum eta for edge masking.
  float edgeMaskMinEta = 1.5;
  /// Threshold for hit sharing between seeds.
  float hitShareThreshold = 0.49;

  // GbtsDataStorage options
  /// Maximum endcap cluster width.
  float maxEndcapClusterWidth = 0.35;
};

}  // namespace Acts::Experimental
