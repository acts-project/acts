// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>

namespace ActsExamples {

/// Print track parameters.
class TrackParametersPrinter : public BareAlgorithm {
 public:
  struct Config {
    /// Input tracks parameters collection.
    std::string inputTrackParameters;
  };

  TrackParametersPrinter(const Config& cfg, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
