// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace FW {

/// Example algorithm that reads/writes data from/to the event store.
class HelloWhiteBoardAlgorithm : public FW::BareAlgorithm {
 public:
  struct Config {
    /// Input collection name.
    std::string input;
    /// Output collection name.
    std::string output;
  };

  HelloWhiteBoardAlgorithm(const Config& cfg,
                           Acts::Logging::Level level = Acts::Logging::INFO);

  /// Read input and copy to the output
  FW::ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace FW
