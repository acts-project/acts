// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/alpaka/utils/queue.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

namespace traccc::alpaka {

/// Helper class for instantiating the correct vecmem objects for the way
/// the @c traccc::alpaka library was built.
class vecmem_objects {
 public:
  /// Constructor from a queue
  explicit vecmem_objects(queue& q);
  /// Move constructor
  vecmem_objects(vecmem_objects&&) noexcept;
  /// Destructor
  ~vecmem_objects();

  /// Move assignment
  vecmem_objects& operator=(vecmem_objects&&) noexcept;

  /// The host memory resource to use
  vecmem::memory_resource& host_mr() const;
  /// The device memory resource to use
  vecmem::memory_resource& device_mr() const;
  /// The shared/managed memory resource to use
  vecmem::memory_resource& shared_mr() const;

  /// The (synchronous) copy object to use
  vecmem::copy& copy() const;
  /// The asynchronous copy object to use
  vecmem::copy& async_copy() const;

 private:
  /// Type holing the implementation
  struct impl;
  /// Smart pointer to the implementation
  std::unique_ptr<impl> m_impl;

};  // class vecmem_objects

}  // namespace traccc::alpaka
