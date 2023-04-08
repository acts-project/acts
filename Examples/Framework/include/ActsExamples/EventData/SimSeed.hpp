// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Seed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <map>
#include <vector>

namespace ActsExamples {
using SimSeed = Acts::Seed<SimSpacePoint>;
/// Container of sim seed
using SimSeedContainer = std::vector<SimSeed>;

}  // namespace ActsExamples
