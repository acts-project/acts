// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/sycl/fitting/kalman_fitting_algorithm.hpp"

namespace traccc::sycl {

kalman_fitting_algorithm::kalman_fitting_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, queue_wrapper& q,
    std::unique_ptr<const Logger> logger, await_function_type await_func)
    : device::kalman_fitting_algorithm{config, mr, copy, std::move(logger)},
      sycl::algorithm_base{q, std::move(await_func)} {}

}  // namespace traccc::sycl
