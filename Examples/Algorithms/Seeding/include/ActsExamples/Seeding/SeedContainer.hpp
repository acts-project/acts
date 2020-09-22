// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Seeding/Seed.hpp"
#include "ActsExamples/Seeding/SimSpacePoint.hpp"

#include <vector>
namespace ActsExamples {
using SeedContainer = std::vector<Acts::Seed<SimSpacePoint> >;
}  // namespace ActsExamples
