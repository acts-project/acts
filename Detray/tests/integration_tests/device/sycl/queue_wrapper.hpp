// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace detray::sycl {

/// Wrapper class for @c sycl::queue
///
/// It is necessary for passing around SYCL queue objects in code that should
/// not be directly exposed to the SYCL headers.
///
/// Note that unlike @c vecmem::sycl::queue_wrapper, this type can not own a
/// queue of its own. It can only view a queue that is owned by "somebody else".
///
class queue_wrapper {
 public:
  /// Wrap an existing @c sycl::queue object, without taking ownership
  queue_wrapper(void* queue) : m_queue(queue) {}

  /// Copy constructor
  queue_wrapper(const queue_wrapper& parent) = default;
  /// Move constructor
  queue_wrapper(queue_wrapper&& parent) = default;

  /// Copy assignment
  queue_wrapper& operator=(const queue_wrapper& rhs) = default;
  /// Move assignment
  queue_wrapper& operator=(queue_wrapper&& rhs) = default;

  /// Access a typeless pointer to the managed @c sycl::queue object
  void* queue() { return m_queue; }
  /// Access a typeless pointer to the managed @c sycl::queue object
  const void* queue() const { return m_queue; }

 private:
  /// Bare pointer to the wrapped @c sycl::queue object
  void* m_queue;

};  // class queue_wrapper

}  // namespace detray::sycl
