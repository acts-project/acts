// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"

namespace ActsExamples {

/// Utility that converts hit coordinates.
class DigitizationCoordinatesConverter final {
 public:
  /// Construct the converter
  ///
  /// @param config is the configuration
  explicit DigitizationCoordinatesConverter(
      DigitizationAlgorithm::Config config);

  /// Get const access to the config
  const DigitizationAlgorithm::Config& config() const { return m_cfg; }

  /// Convert the hit coordinates to the local frame.
  std::tuple<long double, long double> globalToLocal(std::uint64_t moduleId,
                                                     long double x,
                                                     long double y,
                                                     long double z) const;

  /// Convert the hit coordinates to the global frame.
  std::tuple<long double, long double, long double> localToGlobal(
      std::uint64_t moduleId, long double x, long double y) const;

 private:
  /// Configuration
  DigitizationAlgorithm::Config m_cfg;
};

}  // namespace ActsExamples
