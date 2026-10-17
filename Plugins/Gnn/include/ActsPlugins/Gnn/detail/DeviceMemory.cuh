// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Gnn/Tensor.hpp"
#include "ActsPlugins/Gnn/detail/CudaUtils.hpp"

#include <cstddef>
#include <memory_resource>
#include <stdexcept>

#include <cuda_runtime_api.h>
#include <thrust/execution_policy.h>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

namespace ActsPlugins::detail {

/// Non-owning view of the stream, memory resource and copy object of a CUDA
/// execution context, with the accessors of the traccc device algorithms.
/// Memory is released while work on the stream may still use it, so the
/// resource must be stream-ordered (or cache per stream).
class DeviceMemory {
 public:
  /// @throws std::invalid_argument if @p ctx has no memory resource or copy
  explicit DeviceMemory(const ExecutionContext &ctx)
      : m_stream(ctx.stream.value_or(cudaStreamLegacy)),
        m_mr(ctx.memoryResource),
        m_copy(ctx.copy) {
    if (m_mr == nullptr || m_copy == nullptr) {
      throw std::invalid_argument(
          "The CUDA execution context needs a memory resource and a copy "
          "object");
    }
  }

  cudaStream_t stream() const { return m_stream; }
  vecmem::memory_resource &mr() const { return *m_mr; }
  const vecmem::copy &copy() const { return *m_copy; }

  /// Non-synchronizing thrust policy on the stream, temporaries from mr()
  auto policy() const {
    return thrust::cuda::par_nosync(
               std::pmr::polymorphic_allocator<std::byte>(m_mr))
        .on(m_stream);
  }

  void synchronize() const { ACTS_CUDA_CHECK(cudaStreamSynchronize(m_stream)); }

 private:
  cudaStream_t m_stream;
  vecmem::memory_resource *m_mr;
  const vecmem::copy *m_copy;
};

}  // namespace ActsPlugins::detail
