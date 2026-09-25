/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

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
template <typename stream_t, typename value_t = std::byte>
class stream_ordered_allocator {
 public:
  using value_type = value_t;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  stream_ordered_allocator(std::pmr::memory_resource& resource,
                           stream_t& stream) noexcept
      : m_resource(&resource), m_stream(std::addressof(stream)) {}

  template <typename other_value_t>
  stream_ordered_allocator(
      const stream_ordered_allocator<stream_t, other_value_t>& other) noexcept
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
  bool operator==(const stream_ordered_allocator<stream_t, other_value_t>&
                      other) const noexcept {
    return m_stream == other.stream() && *m_resource == *other.resource();
  }

 private:
  std::pmr::memory_resource* m_resource;
  stream_t* m_stream;
};

}  // namespace traccc
