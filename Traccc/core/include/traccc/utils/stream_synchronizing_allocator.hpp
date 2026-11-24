// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <memory_resource>
#include <new>

namespace traccc {

/// Host allocator that waits for stream work before recycling device memory.
///
/// Stream may be any wrapper with a synchronize() member. The wrapper and the
/// resource must outlive the allocator and all its copies/rebound instances.
/// Use the same stream for the execution policy and this allocator.
/// The upstream resource must itself support concurrent host access if shared
/// between threads; stream synchronization only prevents premature reuse.
template <typename value_t, typename stream_t>
class stream_synchronizing_allocator {
 public:
  using value_type = value_t;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  stream_synchronizing_allocator(std::pmr::memory_resource& resource,
                                 stream_t& stream) noexcept
      : m_resource(&resource), m_stream(std::addressof(stream)) {}

  template <typename other_value_t>
  stream_synchronizing_allocator(
      const stream_synchronizing_allocator<other_value_t, stream_t>&
          other) noexcept
      : m_resource(other.resource()), m_stream(other.stream()) {}

  [[nodiscard]] value_type* allocate(size_type n) {
    return static_cast<value_type*>(
        m_resource->allocate(n * sizeof(value_type), alignof(value_type)));
  }

  void deallocate(value_type* ptr, size_type n) {
    // Do not return the allocation to the cache if synchronization fails.
    m_stream->synchronize();
    m_resource->deallocate(ptr, n * sizeof(value_type), alignof(value_type));
  }

  std::pmr::memory_resource* resource() const noexcept { return m_resource; }
  stream_t* stream() const noexcept { return m_stream; }

  template <typename other_value_t>
  bool operator==(const stream_synchronizing_allocator<other_value_t, stream_t>&
                      other) const noexcept {
    return m_stream == other.stream() && *m_resource == *other.resource();
  }

 private:
  std::pmr::memory_resource* m_resource;
  stream_t* m_stream;
};

template <typename stream_t>
stream_synchronizing_allocator(std::pmr::memory_resource&, stream_t&)
    -> stream_synchronizing_allocator<std::byte, stream_t>;

}  // namespace traccc
