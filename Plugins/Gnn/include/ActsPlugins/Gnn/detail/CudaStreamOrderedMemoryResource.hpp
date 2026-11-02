// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Gnn/detail/CudaUtils.hpp"

#include <cstddef>

#include <cuda_runtime_api.h>
#include <vecmem/memory/memory_resource.hpp>

namespace ActsPlugins::detail {

/// vecmem memory resource on the stream-ordered CUDA allocator
/// (cudaMallocAsync / cudaFreeAsync). Freed memory is only reused by work
/// ordered after the free, and the driver pool is thread safe.
class CudaStreamOrderedMemoryResource final : public vecmem::memory_resource {
 public:
  explicit CudaStreamOrderedMemoryResource(cudaStream_t stream)
      : m_stream(stream) {}

  cudaStream_t stream() const { return m_stream; }

 private:
  void *do_allocate(std::size_t bytes, std::size_t /*alignment*/) override {
    void *ptr = nullptr;
    if (bytes > 0) {
      ACTS_CUDA_CHECK(cudaMallocAsync(&ptr, bytes, m_stream));
    }
    return ptr;
  }

  void do_deallocate(void *ptr, std::size_t /*bytes*/,
                     std::size_t /*alignment*/) override {
    // Must not throw (called from destructors); a failure is reported by the
    // next checked CUDA call
    if (ptr != nullptr) {
      static_cast<void>(cudaFreeAsync(ptr, m_stream));
    }
  }

  bool do_is_equal(
      const vecmem::memory_resource &other) const noexcept override {
    const auto *o =
        dynamic_cast<const CudaStreamOrderedMemoryResource *>(&other);
    return o != nullptr && o->m_stream == m_stream;
  }

  cudaStream_t m_stream;
};

}  // namespace ActsPlugins::detail
