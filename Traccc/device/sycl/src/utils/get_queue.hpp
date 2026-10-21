// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/sycl/utils/queue_wrapper.hpp"

// SYCL include(s).
#include <sycl/sycl.hpp>

namespace traccc::sycl::details {

/// Helper function for getting a @c sycl::queue out of
/// @c traccc::sycl::queue_wrapper (non-const)
::sycl::queue& get_queue(traccc::sycl::queue_wrapper& queue);

/// Helper function for getting a @c sycl::queue out of
/// @c traccc::sycl::queue_wrapper (const)
const ::sycl::queue& get_queue(const traccc::sycl::queue_wrapper& queue);

}  // namespace traccc::sycl::details
