// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {
namespace Options {

void addTruthSeedSelectorOptions(Options::Description& desc);

ActsExamples::TruthSeedSelector::Config readTruthSeedSelectorConfig(
    const Options::Variables& vars);

}  // namespace Options
}  // namespace ActsExamples
