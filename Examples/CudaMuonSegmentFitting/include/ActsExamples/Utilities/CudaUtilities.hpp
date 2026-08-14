// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <cuda_runtime_api.h>

namespace ActsExamples {

inline void cudaAssert(cudaError_t code, const char* file, int line) {
  if (code != cudaSuccess) {
    std::stringstream ss;
    ss << "CUDA error: " << cudaGetErrorString(code) << ", " << file << ":"
       << line;
    throw std::runtime_error(ss.str());
  }
}

}  // namespace ActsExamples

#define ACTS_CUDA_CHECK(ans)                             \
  do {                                                   \
    ActsExamples::cudaAssert((ans), __FILE__, __LINE__); \
  } while (0)

namespace ActsExamples {

/// Owning CUDA stream with non-throwing cleanup.
class CudaStream {
 public:
  CudaStream() {
    ACTS_CUDA_CHECK(cudaStreamCreateWithFlags(&m_stream, cudaStreamNonBlocking));
  }

  CudaStream(const CudaStream&) = delete;
  CudaStream& operator=(const CudaStream&) = delete;

  CudaStream(CudaStream&& other) noexcept
      : m_stream{std::exchange(other.m_stream, nullptr)} {}

  CudaStream& operator=(CudaStream&& other) noexcept {
    if (this != &other) {
      reset();
      m_stream = std::exchange(other.m_stream, nullptr);
    }
    return *this;
  }

  ~CudaStream() noexcept { reset(); }

  cudaStream_t get() const noexcept { return m_stream; }

  void synchronize() const { ACTS_CUDA_CHECK(cudaStreamSynchronize(m_stream)); }

 private:
  void reset() noexcept {
    if (m_stream != nullptr) {
      (void)cudaStreamDestroy(m_stream);
      m_stream = nullptr;
    }
  }

  cudaStream_t m_stream = nullptr;
};

template <typename T>
void allocateDeviceColumn(T*& deviceColumn, std::size_t size) {
  if (size == 0) {
    return;
  }

  ACTS_CUDA_CHECK(
      cudaMalloc(reinterpret_cast<void**>(&deviceColumn), size * sizeof(T)));
}

template <typename T>
void freeDeviceColumn(T*& deviceColumn) noexcept {
  if (deviceColumn != nullptr) {
    // Cleanup is used by noexcept destructors and therefore cannot report a
    // cudaFree failure by throwing.
    (void)cudaFree(deviceColumn);
    deviceColumn = nullptr;
  }
}

template <typename T>
void copyColumnToDevice(T* deviceColumn, const std::vector<T>& hostColumn,
                        cudaStream_t stream) {
  if (hostColumn.empty()) {
    return;
  }

  ACTS_CUDA_CHECK(cudaMemcpyAsync(deviceColumn, hostColumn.data(),
                                  hostColumn.size() * sizeof(T),
                                  cudaMemcpyHostToDevice, stream));
}

template <typename T>
void copyColumnToHost(std::vector<T>& hostColumn, const T* deviceColumn,
                      cudaStream_t stream) {
  if (hostColumn.empty() || deviceColumn == nullptr) {
    return;
  }

  ACTS_CUDA_CHECK(cudaMemcpyAsync(hostColumn.data(), deviceColumn,
                                  hostColumn.size() * sizeof(T),
                                  cudaMemcpyDeviceToHost, stream));
}

}  // namespace ActsExamples
