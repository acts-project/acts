// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/hip/utils/stream_wrapper.hpp"

// Vecmem include(s).
#include <vecmem/utils/abstract_event.hpp>

namespace traccc::hip {

/// Suspend execution of a TBB task until work on the HIP stream has completed.
/// A callback-based mechanism is used to handle the completion notification.
///
/// @param event The event to synchronize, unused.
/// @param stream The HIP stream.
///
/// @note Should be called only from a TBB task.
///
void tbb_await_callback(vecmem::abstract_event& event,
                        const stream_wrapper& stream);
}  // namespace traccc::hip
