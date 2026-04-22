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

// SYCL include(s).
#include <sycl/sycl.hpp>

// System include(s).
#include <cstddef>

namespace detray::test::sycl {

/// Execute a test functor using SYCL, on @c array_sizes threads
template <class functor_t, class... Args>
void execute_sycl_test(::sycl::queue& queue, std::size_t array_sizes,
                       Args... args) {
  // Submit a kernel that would run the specified functor.
  queue
      .submit([&](::sycl::handler& h) {
        // Use parallel_for without specifying a "kernel class" explicitly.
        // Unfortunately the functor_t class is too complicated, and DPC++
        // dies on it. While providing a unique simple class for every
        // template specialisation is also pretty impossible. :-(
        h.parallel_for(::sycl::range<1>(array_sizes), [=](::sycl::item<1> id) {
          // Find the current index that we need to
          // process.
          const std::size_t i = id[0];
          if (i >= array_sizes) {
            return;
          }
          // Execute the test functor for this index. Note
          // that std::forward cannot be used here, as the
          // function arguments are passed through as
          // copies to the SYCL kernel. So perfect
          // forwarding is out of question.
          functor_t()(i, args...);
        });
      })
      .wait_and_throw();
}

}  // namespace detray::test::sycl
