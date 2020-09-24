// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// SYCL plugin include(s)
#include "Acts/Plugins/Sycl/Utilities/CalculateNdRange.hpp"

namespace Acts::Sycl {
cl::sycl::nd_range<1> calculate1DimNDRange(const uint32_t numThreads,
                                           const uint32_t workGroupSize) {
  auto q = (numThreads + workGroupSize - 1) / workGroupSize;
  return cl::sycl::nd_range<1>{cl::sycl::range<1>(q * workGroupSize),
                               cl::sycl::range<1>(workGroupSize)};
}

cl::sycl::nd_range<2> calculate2DimNDRange(const uint32_t numThreadsDim0,
                                           const uint32_t numThreadsDim1,
                                           const uint32_t workGroupSize) {
  uint32_t wgSizeDim0 = 1;
  uint32_t wgSizeDim1 = workGroupSize;

  while (numThreadsDim1 < wgSizeDim1 && 1 < wgSizeDim1 && wgSizeDim1 % 2 == 0) {
    wgSizeDim1 /= 2;
    wgSizeDim0 *= 2;
  }
  auto q1 = (numThreadsDim0 + wgSizeDim0 + 1) / wgSizeDim0;
  auto q2 = (numThreadsDim1 + wgSizeDim1 + 1) / wgSizeDim1;
  return cl::sycl::nd_range<2>{
      cl::sycl::range<2>(q1 * wgSizeDim0, q2 * wgSizeDim1),
      cl::sycl::range<2>(wgSizeDim0, wgSizeDim1)};
}
}  // namespace Acts::Sycl
