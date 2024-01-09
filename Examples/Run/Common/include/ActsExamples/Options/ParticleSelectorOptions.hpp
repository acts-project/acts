// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace Options {

void addParticleSelectorOptions(Options::Description& desc);

ActsExamples::ParticleSelector::Config readParticleSelectorConfig(
    const Options::Variables& vars);

}  // namespace Options
}  // namespace ActsExamples
