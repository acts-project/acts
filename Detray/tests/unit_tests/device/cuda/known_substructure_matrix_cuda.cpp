// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray test include(s)
#include "known_substructure_matrix_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <cstddef>
#include <limits>

using namespace detray;

/// The arithmetic of @c ksm::matrix is covered exhaustively by the host unit
/// test. What cannot be covered there is whether it is callable from device
/// code at all: @c DETRAY_HOST_DEVICE expands to nothing for the host
/// compiler, so a host build says nothing about the annotations. This runs the
/// same computation on both sides and requires the results to be identical.
///
/// The comparison is exact rather than approximate on purpose. Both sides
/// execute the same operations on the same inputs in the same order, so the
/// only way a cell can differ is if the device genuinely did something else --
/// a different code path, or a structural cell that lost its compile-time
/// value. Contraction differences are not a concern for a bitwise-identical
/// instruction sequence over integral-valued doubles.
GTEST_TEST(detray_ksm_cuda, matches_host) {
  constexpr std::size_t n = ksm_cuda_test::result_size;

  vecmem::cuda::managed_memory_resource mng_mr;
  auto* device_result = static_cast<double*>(
      mng_mr.allocate(n * sizeof(double), alignof(double)));
  ASSERT_NE(device_result, nullptr);

  for (std::size_t i = 0u; i < n; ++i) {
    device_result[i] = std::numeric_limits<double>::signaling_NaN();
  }

  // Same definition, compiled by the host compiler
  std::array<double, n> host_result{};
  ksm_cuda_test::compute(host_result.data());

  // ... and by nvcc, inside a kernel
  ksm_run_on_device(device_result);

  for (std::size_t i = 0u; i < n; ++i) {
    EXPECT_EQ(host_result.at(i), device_result[i])
        << "host and device disagree at result index " << i;
  }

  mng_mr.deallocate(device_result, n * sizeof(double), alignof(double));
}
