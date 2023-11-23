// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding/TrigBase.hpp"  //definition of Trigsispacepoint base and trigtriplets

#include <memory>

// core algorithm so in acts namespace
namespace Acts {

template <typename T>
class SeedFilter;

template <typename SpacePoint>
struct SeedFinderFTFConfig {
  // // how many sigmas of scattering angle should be considered?
  float sigmaScattering = 5;

  // Seed cut
  float minPt = 400. * Acts::UnitConstants::MeV;

  ///////////some declared not filled in by reco: //////
  std::shared_ptr<Acts::SeedFilter<SpacePoint>> seedFilter;

  //   //detector ROI
  //   // derived values, set on SeedFinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;
  //   bool isInInternalUnits = false;
  /// for load space points
  unsigned int maxSeedsPerSpM = 5;

  float m_phiSliceWidth{};
  float m_nMaxPhiSlice = 53;  // set in athena TrigCombinatorialsettings
  bool m_useClusterWidth = false;
  std::string fastrack_input_file;
  std::vector<TrigInDetSiLayer> m_layerGeometry;

  // for run function
  // m_settings:
  bool m_LRTmode = true;  // eventually want to set from full chain
  bool m_useEtaBinning = true;
  bool m_doubletFilterRZ = true;
  float m_minDeltaRadius = 2.0;  // athena itk 2.0
  float m_tripletD0Max = 4.0;    // m_settings
  unsigned int m_maxTripletBufferLength = 3;

  // from Run function //used to be const but cant as FTF config needs to be
  // changable
  int MaxEdges = 2000000;
  float cut_dphi_max = 0.012;
  float cut_dcurv_max = 0.001;
  float cut_tau_ratio_max = 0.007;
  float maxOuterRadius = 550.0;
  float m_PtMin = 1000.0;
  float m_tripletPtMinFrac = 0.3;
  float m_tripletPtMin =
      m_PtMin *
      m_tripletPtMinFrac;  // used to be 2500, changing this helps low pt a lot
  double ptCoeff = 0.29997 * 1.9972 / 2.0;  // from athena

  // ROI:
  bool containsPhi() {
    return false;
    // need to implement this function
  }

  ////
  // 2 member functions
  SeedFinderFTFConfig calculateDerivedQuantities() const {
    // thorw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderFTFConfig config = *this;
    // use a formula to calculate scattering

    return config;
  }

  SeedFinderFTFConfig toInternalUnits() const {
    // throw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderFTFConfig config = *this;
    // divides inputs by 1mm, all ones input
    // changes member inInInternalUnits to true
    return config;
  }

};  // end of config struct

}  // namespace Acts
