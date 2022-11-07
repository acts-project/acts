// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <limits>
#include <string>

namespace ActsExamples {

/// Select tracks by applying some selection cuts.
class TrackModifier final : public BareAlgorithm {
 public:
  struct Config {
    /// Input track parameters collection.
    std::string inputTrackParameters;
    /// Output track parameters collection.
    std::string outputTrackParameters;

    /// When turned on, only keed the diagonal of the cov matrix.
    bool dropCovariance{false};
    /// Scale cov matrix;
    double covScale{1};
  };

  TrackModifier(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
