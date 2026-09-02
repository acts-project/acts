// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"

#include <cstdlib>
#include <iostream>
#include <optional>
#include <string_view>

namespace ActsTests {

namespace {

std::optional<Acts::Logging::Level> parseLevel(std::string_view name) {
  using Acts::Logging::Level;
  for (Level level : {Level::VERBOSE, Level::DEBUG, Level::INFO, Level::WARNING,
                      Level::ERROR, Level::FATAL, Level::MAX}) {
    if (Acts::Logging::levelName(level) == name) {
      return level;
    }
  }
  return std::nullopt;
}

/// Arms the whole test executable from ACTS_LOG_FAILURE_THRESHOLD.
///
/// Acts::Core does not read the environment, so the tests do it for
/// themselves. Every unit test links this library and it is shared, so this
/// runs once per executable, ahead of the executable's own initialisers.
struct FailureThresholdFromEnvironment {
  FailureThresholdFromEnvironment() {
    const char* value = std::getenv("ACTS_LOG_FAILURE_THRESHOLD");
    if (value == nullptr) {
      return;
    }
    if (std::optional<Acts::Logging::Level> level = parseLevel(value);
        level.has_value()) {
      Acts::Logging::setFailureThreshold(*level);
    } else {
      std::cerr << "ACTS_LOG_FAILURE_THRESHOLD is set to unknown value: "
                << value << std::endl;
    }
  }
};

const FailureThresholdFromEnvironment s_failureThresholdFromEnvironment;

}  // namespace

}  // namespace ActsTests
