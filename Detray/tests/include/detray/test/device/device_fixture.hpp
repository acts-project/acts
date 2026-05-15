// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/utils/approximately_equal.hpp"

// Test include(s).
#include "detray/test/framework/fixture_base.hpp"
#include "detray/test/framework/types.hpp"

// VecMem include(s).
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cmath>
#include <memory>
#include <type_traits>

namespace detray::test {

/// Data fixture for device test cases
template <typename R>
  requires std::is_default_constructible_v<R>
class device_fixture : public fixture_base<> {
 public:
  /// Constructor, setting up the vecmem resource
  device_fixture(vecmem::memory_resource& mr) : m_resource(mr) {}

 protected:
  /// Set up the host and device results
  virtual void SetUp() override {
    m_output_host =
        std::make_unique<vecmem::vector<R>>(this->size(), &m_resource);
    m_output_device =
        std::make_unique<vecmem::vector<R>>(this->size(), &m_resource);

    for (std::size_t i = 0u; i < this->size(); ++i) {
      m_output_host->at(i) = {};
      m_output_device->at(i) = {};
    }
  }

  /// Function resetting the output after the test
  virtual void TearDown() override {
    m_output_host.reset();
    m_output_device.reset();
  }

  /// @returns the number of results
  constexpr std::size_t size() const { return m_size; }

  /// @returns access to the memory resource of the device test suite
  constexpr vecmem::memory_resource& resource() { return m_resource; }

  /// Compare the outputs, after the data processing is finished.
  void compareOutputs() const {
    for (std::size_t i = 0u; i < this->size(); ++i) {
      // Use the inbuilt comparator for algebra types
      if constexpr (detray::concepts::value<R> || detray::concepts::scalar<R> ||
                    detray::concepts::vector<R> ||
                    detray::concepts::matrix<R>) {
        EXPECT_TRUE(detray::algebra::approx_equal(
            m_output_host->at(i), m_output_device->at(i), this->tolerance()))
            << "error at i = " << i << "\nHOST:\n"
            << m_output_host->at(i) << "\nDEVICE:\n"
            << m_output_device->at(i);

      } else {
        EXPECT_EQ(m_output_host->at(i), m_output_device->at(i))
            << "error at i = " << i << "\nHOST:\n"
            << m_output_host->at(i) << "\nDEVICE:\n"
            << m_output_device->at(i);
      }
    }
  }

  /// Memory resource provided by the derived class
  vecmem::memory_resource& m_resource;

  /// Size for the tested arrays.
  static constexpr std::size_t m_size = 5000u;

  /// @name Outputs for the tests
  /// @{
  std::unique_ptr<vecmem::vector<R>> m_output_host;
  std::unique_ptr<vecmem::vector<R>> m_output_device;
  /// @}

};  // class device_fixture

}  // namespace detray::test
