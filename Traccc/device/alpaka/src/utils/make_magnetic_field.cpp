// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/alpaka/utils/make_magnetic_field.hpp"

#include "magnetic_field_types.hpp"

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
#include "traccc/cuda/utils/make_magnetic_field.hpp"
#elif defined(ALPAKA_ACC_GPU_HIP_ENABLED)
#include "traccc/hip/utils/make_magnetic_field.hpp"
#elif defined(ALPAKA_ACC_SYCL_ENABLED)
#include "../utils/get_queue.hpp"
#include "traccc/sycl/utils/make_magnetic_field.hpp"
#endif

// System include(s).
#include <stdexcept>

namespace traccc::alpaka {

magnetic_field make_magnetic_field(const magnetic_field& bfield,
                                   [[maybe_unused]] const queue& queue) {
#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
  return traccc::cuda::make_magnetic_field(
      bfield, traccc::cuda::magnetic_field_storage::global_memory);
#elif defined(ALPAKA_ACC_GPU_HIP_ENABLED)
  return traccc::hip::make_magnetic_field(
      bfield, traccc::hip::magnetic_field_storage::global_memory);
#elif defined(ALPAKA_ACC_SYCL_ENABLED)
  ::sycl::queue q(::alpaka::getNativeHandle(details::get_queue(queue)));
  traccc::sycl::queue_wrapper qw{&q};
  return traccc::sycl::make_magnetic_field(bfield, qw);
#else
  return bfield;
#endif
}

}  // namespace traccc::alpaka
