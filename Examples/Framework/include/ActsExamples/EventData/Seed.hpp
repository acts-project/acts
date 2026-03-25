// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/Types.hpp"

namespace ActsExamples {

using SeedIndex = Acts::SeedIndex2;

using SeedProxy = Acts::SeedContainer2::MutableProxy;
using ConstSeedProxy = Acts::SeedContainer2::ConstProxy;

using SeedContainer = Acts::SeedContainer2;

}  // namespace ActsExamples
