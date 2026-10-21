// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <sstream>
#include <stdexcept>
#include <utility>

#include <cuda_runtime_api.h>

namespace ActsExamples {

/// Owning CUDA stream with non-throwing cleanup.
class CudaStream {
 public:
  CudaStream() {
    check(cudaStreamCreateWithFlags(&m_stream, cudaStreamNonBlocking), __FILE__,
          __LINE__);
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

  void synchronize() const {
    check(cudaStreamSynchronize(m_stream), __FILE__, __LINE__);
  }

 private:
  static void check(cudaError_t code, const char* file, int line) {
    if (code != cudaSuccess) {
      std::stringstream ss;
      ss << "CUDA error: " << cudaGetErrorString(code) << ", " << file << ":"
         << line;
      throw std::runtime_error(ss.str());
    }
  }

  void reset() noexcept {
    if (m_stream != nullptr) {
      (void)cudaStreamDestroy(m_stream);
      m_stream = nullptr;
    }
  }

  cudaStream_t m_stream = nullptr;
};

}  // namespace ActsExamples
