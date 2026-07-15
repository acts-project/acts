/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/device/global_index.hpp"

// SYCL include(s).
#include <sycl/sycl.hpp>

namespace traccc::sycl::details {

/// Function creating a global index in a 1D SYCL kernel
inline device::global_index_t global_index(const ::sycl::nd_item<1>& item) {

    return static_cast<device::global_index_t>(item.get_global_linear_id());
}

/// Function creating a global index in a 2D SYCL kernel
inline device::global_index_t global_index(const ::sycl::nd_item<2>& item) {

    return static_cast<device::global_index_t>(item.get_global_linear_id());
}

/// Function creating a global index in a 3D SYCL kernel
inline device::global_index_t global_index(const ::sycl::nd_item<3>& item) {

    return static_cast<device::global_index_t>(item.get_global_linear_id());
}

}  // namespace traccc::sycl::details
