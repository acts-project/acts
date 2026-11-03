// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/sycl/finding/combinatorial_kalman_filter_algorithm.hpp"

namespace traccc::sycl {

combinatorial_kalman_filter_algorithm::combinatorial_kalman_filter_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, queue_wrapper& queue,
    std::unique_ptr<const Logger> logger,
    std::unique_ptr<traccc::sycl::kalman_fitting_algorithm> kf_fitter,
    await_function_type await_func)
    : device::combinatorial_kalman_filter_algorithm(
          config, mr, copy, std::move(logger), std::move(kf_fitter)),
      sycl::algorithm_base(queue, std::move(await_func)) {}

}  // namespace traccc::sycl
