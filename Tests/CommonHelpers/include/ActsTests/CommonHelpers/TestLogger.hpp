// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <iostream>
#include <memory>
#include <string>

namespace ActsTests {

/// The level at or above which a test logger fails the test.
inline constexpr Acts::Logging::Level s_testFailureThreshold =
    Acts::Logging::WARNING;

/// @brief Make a logger for a test, armed to fail on warnings and worse
///
/// A test that logs an expected error disarms its logger with
/// @ref Acts::Logger::withoutFailureThreshold.
///
/// @param name name of the logging instance
/// @param level level below which messages are filtered out
/// @param stream stream to print to
/// @return the armed logger
inline std::unique_ptr<const Acts::Logger> getTestLogger(
    const std::string& name, Acts::Logging::Level level = Acts::Logging::INFO,
    std::ostream* stream = &std::cout) {
  return Acts::getDefaultLogger(name, level, stream, s_testFailureThreshold);
}

}  // namespace ActsTests
