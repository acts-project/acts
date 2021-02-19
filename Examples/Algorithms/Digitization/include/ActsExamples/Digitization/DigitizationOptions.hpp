// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {
namespace Options {

/// Add Digitization options.
///
/// @param desc The options description to add options to
void addDigitizationOptions(Description& desc);

/// Read SmearingAlgorithm Config from the options.
///
/// @param variables The variables to read from
Digitization::AlgorithmConfig readSmearingConfig(const Variables& variables);

/// Read Digitization Config from the options.
///
/// @param variables The variables to read from
Digitization::AlgorithmConfig readDigitizationConfig(const Variables& variables);

Digitization::AlgorithmConfig configureDigitization(const Variables &vm);

std::shared_ptr<ActsExamples::IAlgorithm> createDigitizationAlgorithm(Digitization::AlgorithmConfig &cfg, Acts::Logging::Level lvl);

Acts::GeometryHierarchyMap<DigitizationConfig> readConfigFromJson(const std::string &path);

std::vector<std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
getBoundIndices(ActsExamples::Digitization::AlgorithmConfig &cfg);


}  // namespace Options
}  // namespace ActsExamples
