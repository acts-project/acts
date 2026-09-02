// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

namespace ActsExamples {

/// @brief Get the threshold the framework arms the loggers of a job at
///
/// Sequence elements build their own loggers from a level, so a job cannot
/// hand them one. This is the value they are armed at instead.
///
/// @return the threshold, or @ref Acts::Logging::Level::MAX for "do not arm"
Acts::Logging::Level getLogFailureThreshold();

/// @brief Set the threshold the framework arms the loggers of a job at
///
/// @warning Global state, and only reaches loggers built after the call.
/// @param level The new threshold, or @ref Acts::Logging::Level::MAX to disarm
void setLogFailureThreshold(Acts::Logging::Level level);

}  // namespace ActsExamples
