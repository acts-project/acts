// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"

namespace ActsExamples {

using SpacePointIndex = Acts::SpacePointIndex2;

using SpacePointProxy = Acts::SpacePointContainer2::MutableProxy;
using ConstSpacePointProxy = Acts::SpacePointContainer2::ConstProxy;

using SpacePointColumns = Acts::SpacePointColumns;
using SpacePointContainer = Acts::SpacePointContainer2;

}  // namespace ActsExamples
