// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/alpaka/utils/queue.hpp"
#include "utils.hpp"

namespace traccc::alpaka::details {

/// Helper function for getting a @c Queue out of @c queue (non-const)
Queue& get_queue(queue& q);

/// Helper function for getting a @c Queue out of @c queue (const)
const Queue& get_queue(const queue& q);

}  // namespace traccc::alpaka::details
