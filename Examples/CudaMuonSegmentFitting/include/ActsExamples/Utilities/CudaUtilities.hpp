// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <sstream>

#include <cuda_runtime_api.h>

namespace ActsExamples {

inline void cudaAssert(cudaError_t code, const char *file, int line) {
  if (code != cudaSuccess) {
    std::stringstream ss;
    ss << "CUDA error: " << cudaGetErrorString(code) << ", " << file << ":"
       << line;
    throw std::runtime_error(ss.str());
  }
}

}  // namespace ActsExamples

#define ACTS_CUDA_CHECK(ans)                                    \
  do {                                                          \
    ActsExamples::cudaAssert((ans), __FILE__, __LINE__); \
  } while (0)

namespace ActsExamples {

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
    ACTS_CUDA_CHECK(cudaFree(deviceColumn));
    deviceColumn = nullptr;
  }
}

template <typename T>
void copyColumnToDevice(T* deviceColumn, const std::vector<T>& hostColumn) {
  if (hostColumn.empty()) {
    return;
  }

  ACTS_CUDA_CHECK(cudaMemcpy(deviceColumn, hostColumn.data(),
                       hostColumn.size() * sizeof(T), cudaMemcpyHostToDevice));
}

template <typename T>
void copyColumnToHost(std::vector<T>& hostColumn, const T* deviceColumn) {
  if (hostColumn.empty() || deviceColumn == nullptr) {
    return;
  }

  ACTS_CUDA_CHECK(cudaMemcpy(hostColumn.data(), deviceColumn,
                       hostColumn.size() * sizeof(T), cudaMemcpyDeviceToHost));
}

}  // namespace ActsExamples
