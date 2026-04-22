// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <cstddef>
#include <utility>

namespace detray::test {

/// Execute a test functor on the host (later compare with device results)
template <class functor_t, class... Args>
void execute_host_test(std::size_t array_sizes, Args... args) {
  // Instantiate the functor.
  constexpr functor_t functor;

  // Execute the functor on all elements of the array(s).
  for (std::size_t i = 0; i < array_sizes; ++i) {
    functor(i, std::forward<Args>(args)...);
  }
}

}  // namespace detray::test
