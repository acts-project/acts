// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Fatras/FatrasSimulation.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {
namespace Options {

void addFatrasOptions(Options::Description& desc);

ActsExamples::FatrasSimulation::Config readFatrasConfig(
    const Options::Variables& vars);

}  // namespace Options
}  // namespace ActsExamples
