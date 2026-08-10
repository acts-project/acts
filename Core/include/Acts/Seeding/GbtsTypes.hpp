// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"

namespace Acts::Experimental {

/// Index of a graph node inside GbtsNodeStorage. Nodes are stored ordered by
/// (eta bin, phi), so all nodes of an eta bin form a contiguous range.
using GbtsNodeIndex = SpacePointIndex;

/// Sentinel for an unset graph node index.
static constexpr GbtsNodeIndex kGbtsNodeIndexInvalid = kSpacePointIndexInvalid;

}  // namespace Acts::Experimental
