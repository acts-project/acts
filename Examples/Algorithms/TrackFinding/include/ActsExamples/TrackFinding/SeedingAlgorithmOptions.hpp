// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {
namespace Options {

/// Add SeedingAlgorithm options.
///
/// @param desc The options description to add options to
void addSeedingAlgorithmOptions(Description& desc);

//Getting printouts for EA
void addMLOutputOptions(Description& desc);

/// Read SeedingAlgorithm options to create the algorithm config.
///
/// @param variables The variables to read from
SeedingAlgorithm::Config readSeedingAlgorithmConfig(const Variables& variables);

bool readMLOutputConfig(const Variables& variables);

}  // namespace Options
}  // namespace ActsExamples
