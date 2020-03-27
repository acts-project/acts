// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <string>

#include "ACTFW/Framework/BareService.hpp"

namespace FW {

/// A simple service that adds an event block index.
class HelloService : public BareService {
 public:
  struct Config {
    /// How many events in one block.
    std::size_t eventsPerBlock = 4;
    /// Under which name to store the block index.
    std::string blockIndexName = "eventBlock";
  };

  HelloService(const Config& cfg, Acts::Logging::Level level);

  void startRun() final override;

  void prepare(AlgorithmContext& ctx) final override;

 private:
  Config m_cfg;
};

}  // namespace FW
