// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// SYCL include(s).
#include <CL/sycl.hpp>

// System include(s).
#include <cstddef>
#include <memory>

namespace Acts::Sycl {

namespace detail {

/// Deleter functor for the smart pointer type(s)
class DeviceArrayDeleter {
 public:
  /// Constructor, with the queue that the memory was allocated on/with
  DeviceArrayDeleter(cl::sycl::queue& queue) : m_queue(&queue) {}
  /// Operator performing the deletion of the memory
  void operator()(void* ptr) {
    if (ptr != nullptr) {
      cl::sycl::free(ptr, *m_queue);
    }
  }

 private:
  /// The queue that manages the memory area in question
  cl::sycl::queue* m_queue;
};  // class DeviceArrayDeleter

}  // namespace detail

/// Convenience type for using (primitive) variable arrays on a SYCL device
template <typename T>
using device_array = std::unique_ptr<T, detail::DeviceArrayDeleter>;

/// Function creating a primitive array in SYCL device memory
template <typename T>
device_array<T> make_device_array(size_t size, cl::sycl::queue& queue) {
  return device_array<T>(cl::sycl::malloc_device<T>(size, queue),
                         detail::DeviceArrayDeleter(queue));
}

}  // namespace Acts::Sycl
