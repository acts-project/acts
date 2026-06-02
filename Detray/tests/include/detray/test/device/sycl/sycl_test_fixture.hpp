// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Test include(s).
#include "detray/test/device/device_fixture.hpp"

// VecMem include(s).
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/sycl/shared_memory_resource.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

// SYCL include(s).
#include <sycl/sycl.hpp>

namespace detray::test::sycl {

/// SYCL specific test fixture class
template <class base_fixture_t>
class sycl_test_fixture : public base_fixture_t {
 public:
  /// Constructor, setting up the inputs for all of the tests
  sycl_test_fixture() : base_fixture_t(m_resource) {}

 protected:
  /// Queue to be used by all of the tests.
  ::sycl::queue m_queue;
  /// Memory resource for all of the tests.
  vecmem::sycl::shared_memory_resource m_resource{&m_queue};
};

}  // namespace detray::test::sycl
