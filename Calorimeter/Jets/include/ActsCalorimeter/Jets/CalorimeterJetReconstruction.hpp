// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsCalorimeter/EventData/CalorimeterCluster.hpp"
#include "ActsCalorimeter/EventData/CalorimeterJet.hpp"

#include <span>

namespace ActsCalorimeter {

/// Inclusive anti-kt jets from calibrated cluster energies and directions.
/// Configuration is immutable and all clustering state belongs to one call.
class CalorimeterJetReconstruction {
 public:
  struct Config {
    /// FastJet radius in rapidity-phi space.
    double radius = 0.4;
    /// Inclusive minimum jet transverse momentum in ACTS native units.
    double jetPtMin = 0;
    /// Fixed reference point for cluster directions, in global ACTS
    /// coordinates. An event-dependent vertex requires a separate adapter in a
    /// later stage.
    Acts::Vector3 origin = Acts::Vector3::Zero();
  };

  explicit CalorimeterJetReconstruction(const Config& config);

  /// Assign p = E * unit(position - origin) to each positive-energy cluster,
  /// then recombine by four-vector addition. Zero-energy clusters are skipped.
  /// Outputs are sorted by descending pT, with constituent indices breaking
  /// exact ties. Constituents are copied before FastJet's sequence is
  /// destroyed. Invalid kinematics or undefined positive-energy directions are
  /// errors.
  CalorimeterJetContainer operator()(
      std::span<const CalorimeterCluster> clusters) const;

 private:
  Config m_config;
};

}  // namespace ActsCalorimeter
