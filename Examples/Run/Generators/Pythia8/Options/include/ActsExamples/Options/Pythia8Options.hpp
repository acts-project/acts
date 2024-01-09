// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {
namespace Options {

/// Add options for a Pythia8-based event-generator.
void addPythia8Options(Description& desc);

/// Create the event generator config from Pythia8 options.
///
/// This builds a full event generator with separate hard scatter and pileup.
/// Not just the Pythia8 process generators to simplify the handling.
EventGenerator::Config readPythia8Options(const Variables& vars,
                                          Acts::Logging::Level lvl);

}  // namespace Options
}  // namespace ActsExamples
