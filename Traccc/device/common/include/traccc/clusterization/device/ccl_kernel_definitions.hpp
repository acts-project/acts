// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// System include(s).
#include <cstddef>

namespace traccc::device::details {
/// These indices in clusterization will only range from 0 to
/// max_cells_per_partition, so we only need a short. Keep a fallback
/// unsigned int index for enormous partitions.
using index_t = unsigned short;
using fallback_index_t = unsigned int;

/// The limit on the stack size in terms of cells per thread.
static constexpr std::size_t CELLS_PER_THREAD_STACK_LIMIT = 32;
}  // namespace traccc::device::details
