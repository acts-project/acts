// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "traccc/alpaka/utils/await.hpp"

// System include(s).
#include <concepts>

namespace traccc::alpaka {

void await_sync_event(vecmem::abstract_event& event, const queue& /*queue*/) {
  event.wait();
}

static_assert(
    std::constructible_from<await_function_type, decltype(await_sync_event)>,
    "await_sync_event should be of type await_function_type");

}  // namespace traccc::alpaka
