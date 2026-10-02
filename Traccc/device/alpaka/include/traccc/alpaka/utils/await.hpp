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

// VecMem include(s).
#include <vecmem/utils/abstract_event.hpp>

// System include(s).
#include <functional>

namespace traccc::alpaka {

/// Type of the function used to synchronize async operations.
using await_function_type =
    std::function<void(vecmem::abstract_event&, const queue&)>;

/// Synchronize an event by waiting for it to complete.
///
/// @param event The event to synchronize.
/// @param queue The Alpaka queue, unused.
///
void await_sync_event(vecmem::abstract_event& event, const queue& queue);

}  // namespace traccc::alpaka
