// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/Generators/EventGenerator.hpp"
#include "ACTFW/Utilities/OptionsFwd.hpp"

namespace FW {
namespace Options {

/// Options for a Pythia8-based event-generator w/ hard scatter, variable
/// pile-up, and smeared vertices.
void addPythia8Options(boost::program_options::options_description& opt);

/// Create the event generator config from the options
///
/// This builds a full event generator with separate hard scatter and pileup.
/// Not just the Pythia8 process generators to simplify the handling.
EventGenerator::Config readPythia8Options(
    const boost::program_options::variables_map& vm, Acts::Logging::Level lvl);

}  // namespace Options
}  // namespace FW
