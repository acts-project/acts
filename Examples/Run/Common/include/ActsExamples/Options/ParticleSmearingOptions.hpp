// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <boost/program_options.hpp>

namespace ActsExamples {

class Sequencer;

namespace Options {

/// Add particle smearing options with a `smearing-` prefix.
void addParticleSmearingOptions(Description& desc);

/// Read options
ParticleSmearing::Config readParticleSmearingOptions(
    const Options::Variables& vars);

}  // namespace Options
}  // namespace ActsExamples
