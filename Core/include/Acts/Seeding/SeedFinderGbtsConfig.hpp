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

/// core algorithm so in acts namespace
namespace Acts::Experimental {

struct SeedFinderGbtsConfig {
  // GbtsSeedingAlgorithm options
  bool BeamSpotCorrection = false;

  // Path to the connector configuration file that defines the layer connections
  std::string connectorInputFile;

  std::string lutInputFile;

  /// SeedFinderGbts option
  bool LRTmode = false;
  bool useML = false;  // use cluster width
  bool matchBeforeCreate = false;
  bool useOldTunings = false;
  float tau_ratio_cut = 0.007;
  float tau_ratio_precut = 0.009f;
  // specify non-zero to override eta bin width from connection file
  // (default 0.2 in createLinkingScheme.py)
  float etaBinOverride = 0.0f;
  // used to calculate phi slices
  float nMaxPhiSlice = 53;
  float minPt = 1.0f * UnitConstants::GeV;
  // derived in CreatSeeds function
  float phiSliceWidth{};

  /// BuildTheGraph() options

  // ~0.3*B/2 - assumes nominal field of 2*T
  double ptCoeff =
      0.29997 * 1.9972 / 2.0 * (UnitConstants::MeV / UnitConstants::mm);
  // bool to use eta binning from geometry structure
  bool useEtaBinning = true;
  // bool applies new Z cuts on doublets
  bool doubletFilterRZ = true;
  // max number of Gbts edges/doublets
  std::int32_t nMaxEdges = 2000000;
  float minDeltaRadius = 2.0;

  /// GbtsTrackingFilter
  // Update()
  float sigmaMS = 0.016;  // for 900 MeV track at eta=0
  float radLen = 0.025;   // 2.5% per layer

  float sigma_x = 0.08;
  float sigma_y = 0.25;

  float weight_x = 0.5;
  float weight_y = 0.5;

  float maxDChi2_x = 5.0;
  float maxDChi2_y = 6.0;

  float add_hit = 14.0;

  float max_curvature = 1e-3f;
  float max_z0 = 170.0;

  // extractSeedsFromTheGraph()
  float edge_mask_min_eta = 1.5;
  float hit_share_threshold = 0.49;

  // GbtsDataStorage
  float max_endcap_clusterwidth = 0.35;

  // 2 member functions
  SeedFinderGbtsConfig calculateDerivedQuantities() const {
    // thorw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderGbtsConfig config = *this;
    // use a formula to calculate scattering

    return config;
  }

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
