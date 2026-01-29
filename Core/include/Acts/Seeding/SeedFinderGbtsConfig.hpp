// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style

#include "Acts/Definitions/Units.hpp"

#include <string>

namespace Acts::Experimental {

/// Configuration options for the GBTs seed finder.
struct SeedFinderGbtsConfig {
  // GbtsSeedingAlgorithm options
  /// Enable beam spot correction.
  bool BeamSpotCorrection = false;

  // Path to the connector configuration file that defines the layer connections
  /// Connector configuration file path.
  std::string connectorInputFile;

  /// Look-up table input file path.
  std::string lutInputFile;

  // SeedFinderGbts option
  /// Enable Large Radius Tracking mode.
  bool lrtMode = false;
  /// Use machine learning features (e.g., cluster width).
  bool useML = false;  // use cluster width
  /// Match seeds before creating them.
  bool matchBeforeCreate = false;
  /// Use legacy tuning parameters.
  bool useOldTunings = false;
  /// Tau ratio cut threshold.
  float tau_ratio_cut = 0.007;
  /// Tau ratio precut threshold.
  float tau_ratio_precut = 0.009f;
  /// Eta bin width override (0 uses default from connection file).
  /// specify non-zero to override eta bin width from connection file
  /// (default 0.2 in createLinkingScheme.py)
  /// Maximum number of phi slices.
  float etaBinOverride = 0.0f;
  /// used to calculate phi slices
  /// Minimum transverse momentum.
  float nMaxPhiSlice = 53;
  float minPt = 1.0f * UnitConstants::GeV;
  /// Phi slice width (derived in CreateSeeds function).
  /// derived in CreatSeeds function
  float phiSliceWidth{};

  // graph building options
  /// Transverse momentum coefficient (~0.3*B/2 - assumes nominal field of 2*T).

  // ~0.3*B/2 - assumes nominal field of 2*T
  double ptCoeff =
      0.29997 * 1.9972 / 2.0;  // ~0.3*B/2 - assumes nominal field of 2*T
  //// Use eta binning from geometry structure.
  bool useEtaBinning = true;  /// bool to use eta binning from geometry structure
  /// Apply RZ cuts on doublets.
  bool doubletFilterRZ = true;  // bool applies new Z cuts on doublets
  //// Maximum number of Gbts edges/doublets.
  std::int32_t nMaxEdges = 2000000;  // max number of Gbts edges/doublets
  /// Minimum delta radius between layers.
  float minDeltaRadius = 2.0;

  // GbtsTrackingFilter options

  /// Multiple scattering sigma (for 900 MeV track at eta=0).
  float sigmaMS = 0.016;  // for 900 MeV track at eta=0
  /// Radiation length fraction per layer (2.5% per layer).
  float radLen = 0.025;  // 2.5% per layer

  /// Measurement uncertainty in x direction.
  float sigma_x = 0.08;
  /// Measurement uncertainty in y direction.
  float sigma_y = 0.25;

  /// Measurement weight in x direction.
  float weight_x = 0.5;
  /// Measurement weight in y direction.
  float weight_y = 0.5;

  /// Maximum delta chi2 in x direction.
  float maxDChi2_x = 5.0;
  /// Maximum delta chi2 in y direction.
  float maxDChi2_y = 6.0;

  /// Chi2 penalty for adding a hit.
  float add_hit = 14.0;

  /// Maximum track curvature.
  float max_curvature = 1e-3f;
  /// Maximum longitudinal impact parameter.
  float max_z0 = 170.0;

  // extractSeedsFromTheGraph()
  /// Minimum eta for edge masking.
  float edge_mask_min_eta = 1.5;
  /// Threshold for hit sharing between seeds.
  float hit_share_threshold = 0.49;

  // GbtsDataStorage
  /// Maximum endcap cluster width.
  float max_endcap_clusterwidth = 0.35;

  // 2 member functions
  /// Calculate derived configuration quantities.
  /// @return Configuration with derived quantities calculated
  SeedFinderGbtsConfig calculateDerivedQuantities() const {
    // thorw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderGbtsConfig config = *this;
    // use a formula to calculate scattering

    return config;
  }

  /// Convert configuration to internal units.
  /// @return Configuration with values in internal units
  SeedFinderGbtsConfig toInternalUnits() const {
    // throw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderGbtsConfig config = *this;
    // divides inputs by 1mm, all ones input
    // changes member inInInternalUnits to true
    return config;
  }

};  // end of config struct

}  // namespace Acts::Experimental
