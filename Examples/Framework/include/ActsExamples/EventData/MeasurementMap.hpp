// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Index.hpp"

#include <unordered_set>
#include <vector>

namespace ActsExamples {

/// Set of measurement indices (positions in MeasurementContainer) that have
/// been consumed by a tracking pass.
///
/// Indices are always relative to the original, unfiltered MeasurementContainer
/// produced by digitization so they remain stable across multiple passes.
using MeasurementMap = std::unordered_set<Index>;

/// Maps a filtered-container index to its index in the original
/// MeasurementContainer.  Element i of this vector holds the original index
/// of the i-th measurement in the filtered container.
using MeasurementIndexRemapping = std::vector<Index>;

}  // namespace ActsExamples
