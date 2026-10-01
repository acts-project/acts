// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <filesystem>
#include <vector>

namespace ActsPlugins {

struct MillePedeParameterResult {
  int label = 0;     ///> label of the parameter
  double val = 0;    ///> fitted value
  double start = 0;  ///> starting value
  double delta = 0;  ///> parameter change in this iteration
  double sigma = 0;  ///> uncertainty
  int nRecords = 0;  ///> number of measurements affected by the parameter
};

/// @brief Read the alignment parameters from a result file.
/// @param mpFile: Location of the file to parse
/// @param logger: logger to print to
Acts::Result<std::vector<MillePedeParameterResult>> readMillePedeResult(
    const std::filesystem::path& mpFile, const Acts::Logger& logger);

}  // namespace ActsPlugins
