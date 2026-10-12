// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsCalorimeter/EventData/CalorimeterHit.hpp"

#include <limits>
#include <span>

namespace ActsCalorimeter {

/// Deterministic prototype response, without noise, smearing or saturation.
class CalorimeterResponse {
 public:
  struct Config {
    /// Dimensionless factor applied once to the sum of accepted deposits.
    double energyScale = 1;
    /// Inclusive threshold in calibrated ACTS energy units.
    double energyThreshold = 0;
    /// Inclusive deposit-time acceptance in ACTS time units.
    double timeMin = -std::numeric_limits<double>::infinity();
    double timeMax = std::numeric_limits<double>::infinity();
  };

  /// Reject invalid scale, threshold or time bounds at construction.
  explicit CalorimeterResponse(const Config& config);

  /// Validate input, select deposits by time, merge by cell, calibrate and
  /// threshold. Output is sorted by cell ID. Zero-energy deposits are ignored.
  /// Accepted deposits with repeated IDs must have identical cell centres.
  /// @throws std::invalid_argument for malformed deposits/cell centres.
  /// @throws std::overflow_error if the accumulated output is not finite.
  CalorimeterHitContainer operator()(
      std::span<const SimCalorimeterHit> deposits) const;

  const Config& config() const { return m_config; }

 private:
  Config m_config;
};

}  // namespace ActsCalorimeter
