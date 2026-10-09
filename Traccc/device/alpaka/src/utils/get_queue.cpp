// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "get_queue.hpp"

// System include(s).
#include <cassert>

namespace traccc::alpaka::details {

Queue& get_queue(queue& q) {
  assert(q.alpakaQueue() != nullptr);
  return *(reinterpret_cast<Queue*>(q.alpakaQueue()));
}

const Queue& get_queue(const queue& q) {
  assert(q.alpakaQueue() != nullptr);
  return *(reinterpret_cast<const Queue*>(q.alpakaQueue()));
}

}  // namespace traccc::alpaka::details
