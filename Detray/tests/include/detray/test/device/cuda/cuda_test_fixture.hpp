// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// VecMem include(s).
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

namespace detray::test::cuda {

/// CUDA specific test fixture class
template <class base_fixture_t>
class cuda_test_fixture : public base_fixture_t {
 public:
  /// Constructor, providing CUDA specific infrastructure to the test fixture
  cuda_test_fixture() : base_fixture_t(m_resource) {}

 protected:
  /// Memory resource for all of the tests.
  vecmem::cuda::managed_memory_resource m_resource;
};

}  // namespace detray::test::cuda
