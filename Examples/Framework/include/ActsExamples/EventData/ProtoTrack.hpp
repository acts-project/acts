// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Index.hpp"

#include <vector>

namespace ActsExamples {

/// A proto track is a collection of hits identified by their indices.
using ProtoTrack = std::vector<Index>;
/// Container of proto tracks. Each proto track is identified by its index.
using ProtoTrackContainer = std::vector<ProtoTrack>;

}  // namespace ActsExamples
