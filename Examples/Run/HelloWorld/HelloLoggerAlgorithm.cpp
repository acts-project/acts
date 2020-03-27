// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "HelloLoggerAlgorithm.hpp"

#include <cstddef>

#include "ACTFW/Framework/WhiteBoard.hpp"

FW::HelloLoggerAlgorithm::HelloLoggerAlgorithm(Acts::Logging::Level level)
    : FW::BareAlgorithm("HelloLogger", level) {}

FW::ProcessCode FW::HelloLoggerAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // using hard-coded data name should be avoided, but i'm a bit lazy tonight.
  auto block = ctx.eventStore.get<std::size_t>("eventBlock");

  ACTS_INFO(" Hello World! (from event=" << ctx.eventNumber
                                         << ", block=" << block << ")");
  ACTS_DEBUG("  - that's an ACTS_DEBUG message");
  ACTS_VERBOSE("  - that's an ACTS_VERBOSE message");
  return FW::ProcessCode::SUCCESS;
}
