// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/sycl/utils/queue_wrapper.hpp"

// VecMem include(s).
#include <vecmem/utils/abstract_event.hpp>

// System include(s).
#include <functional>

namespace traccc::sycl {

/// Type of the function used to synchronize async operations.
using await_function_type =
    std::function<void(vecmem::abstract_event&, const queue_wrapper&)>;

/// Synchronize an event by waiting for it to complete.
///
/// @param event The event to synchronize.
/// @param queue The SYCL queue, unused.
///
void await_sync_event(vecmem::abstract_event& event,
                      const queue_wrapper& /*queue*/);

}  // namespace traccc::sycl
