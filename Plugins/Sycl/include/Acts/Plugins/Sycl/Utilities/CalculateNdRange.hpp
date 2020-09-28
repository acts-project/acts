// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// SYCL include
#include <CL/sycl.hpp>

namespace Acts::Sycl {
/// @brief Calculate global range of 1 dimensional nd_range for
/// kernel execution
///
/// @param[in] numThreads is the number of threads globally
/// @param[in] workGroupSize is the number of threads in one work group
///
/// Calculates the global dimension of the nd_range, which is the smallest
/// multiple of workGroupSize that is not smaller than numThreads.
///
/// @return a one dimensional nd_range of threads
cl::sycl::nd_range<1> calculate1DimNDRange(const uint32_t numThreads,
                                           const uint32_t workGroupSize);

/// @brief Calculate global and local range of 2 dimensional nd_range for
/// kernel execution
///
/// @param[in] numThreadsDim0 is the number of threads globally in the first
/// dimension
/// @param[in] numThreadsDim1 is the number of threads globally in the
/// second dimension
/// @param[in] workGroupSize is the number of threads in one work group
///
/// Local range is calculated the following way:
///  - local range dimensions multiplied together should be equal to
///  workGroupSize
///  - if workGroupSize > numThreadsDim1: set local range in the second
///  dimension to be equal to the smallest factor of 2, that is not smaller
///  that the number of threads globally in the second dimension
///  - else: local range is {1, workGroupSize}
///
/// Global range is calculated the following way:
/// - set the number of threads in both dimensions to the smallest multiple
/// of the work group size in that dimension
///
/// @return a two dimensional nd_range of threads
cl::sycl::nd_range<2> calculate2DimNDRange(const uint32_t numThreadsDim0,
                                           const uint32_t numThreadsDim1,
                                           const uint32_t workGroupSize);

}  // namespace Acts::Sycl
