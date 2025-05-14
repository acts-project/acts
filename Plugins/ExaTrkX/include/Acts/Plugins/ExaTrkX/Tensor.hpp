// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <functional>
#include <numeric>
#include <optional>
#include <ostream>

#ifdef ACTS_EXATRKX_WITH_CUDA
#include <cuda_runtime_api.h>
#else
using cudaStream_t = void *;
#endif

namespace Acts {

/// A simple device description struct
struct Device {
  enum class Type { eCPU, eCUDA };
  Type type = Type::eCPU;
  std::size_t index = 0;

  static Device Cpu() { return {Type::eCPU, 0}; }
  static Device Cuda(std::size_t index = 0) { return {Type::eCUDA, index}; }

  bool isCpu() const { return type == Type::eCPU; }
  bool isCuda() const { return type == Type::eCUDA; }

  bool operator==(const Device &) const = default;
};

/// Easy printout of device
inline std::ostream &operator<<(std::ostream &os, Device device) {
  if (device.type == Device::Type::eCPU) {
    os << "CPU";
  } else {
    os << "CUDA(" << device.index << ")";
  }
  return os;
}

/// Capture the context of the execution
struct ExecutionContext {
  Acts::Device device{Acts::Device::Type::eCPU};
  std::optional<cudaStream_t> stream;
};

/// This is a very small, limited class that models a 2D tensor of arbitrary
/// type. It is move-only, and only possible to create via static factory
/// functions to ensure lifetime management.
/// This on purpose does not implement operations such as clone/to-host/to-cuda
template <typename T>
class Tensor {
 public:
  using Shape = std::array<std::size_t, 2>;
  using Deleter = std::function<void(T *)>;

  static Tensor Create(Shape shape, const ExecutionContext &execContext) {
    T *ptr{};
    Deleter del;

    if (execContext.device.type == Acts::Device::Type::eCPU) {
      ptr = new T[shape[0] * shape[1]];
      del = [](T *p) { delete[] p; };
    } else {
#ifdef ACTS_EXATRKX_WITH_CUDA
      assert(execContext.stream.has_value());
      auto stream = *execContext.stream;
      ACTS_CUDA_CHECK(cudaMallocAsync((void **)&ptr,
                                      sizeof(T) * shape[0] * shape[1], stream));
      del = [stream](T *p) { ACTS_CUDA_CHECK(cudaFreeAsync(p, stream)); };
#else
      throw std::runtime_error(
          "Cannot create CUDA tensor, library was not compiled with CUDA");
#endif
    }
    return Tensor(ptr, shape, del, execContext.device);
  }

  /// Clone the tensor, copying the data to the new device
  /// @param to The {device, stream} to clone to
  /// @note This is a always a deep copy, even if the source and destination are the
  /// same device
  Tensor clone(const ExecutionContext &to) const {
    auto clone = Create(m_shape, to);

    if (m_device.isCpu() && to.device.isCpu()) {
      std::memcpy(clone.data(), m_data, nbytes());
    } else {
#ifdef ACTS_EXATRKX_WITH_CUDA
      assert(to.stream.has_value());
      if (m_device.isCuda() && to.device.isCuda()) {
        ACTS_CUDA_CHECK(cudaMemcpyAsync(clone.data(), m_data, nbytes(),
                                        cudaMemcpyDeviceToDevice, *to.stream));
      } else if (m_device.isCpu() && to.device.isCuda()) {
        ACTS_CUDA_CHECK(cudaMemcpyAsync(clone.data(), m_data, nbytes(),
                                        cudaMemcpyHostToDevice, *to.stream));
      } else if (m_device.isCuda() && to.device.isCpu()) {
        ACTS_CUDA_CHECK(cudaMemcpyAsync(clone.data(), m_data, nbytes(),
                                        cudaMemcpyDeviceToHost, *to.stream));
      }
#else
      throw std::runtime_error(
          "Cannot clone CUDA tensor, library was not compiled with CUDA");
#endif
    }
    return clone;
  }

  Tensor(const Tensor &) = delete;
  Tensor &operator=(const Tensor &) = delete;

  Tensor(Tensor &&other) { moveConstruct(std::move(other)); }

  Tensor &operator=(Tensor &&other) {
    moveConstruct(std::move(other));
    return *this;
  }

  ~Tensor() {
    if (m_deleter) {
      m_deleter(m_data);
    }
  }

  T *data() { return m_data; }
  const T *data() const { return m_data; }
  Shape shape() const { return m_shape; }
  std::size_t size() const { return m_shape[0] * m_shape[1]; }
  std::size_t nbytes() const { return size() * sizeof(T); }
  Acts::Device device() const { return m_device; }

 private:
  Tensor(T *ptr, Shape shape, Deleter deleter, Acts::Device device)
      : m_data(ptr),
        m_shape(shape),
        m_deleter(std::move(deleter)),
        m_device(device) {}

  void moveConstruct(Tensor &&other) {
    // Swap the deleters, so there is no double free
    std::swap(m_deleter, other.m_deleter);
    m_data = other.m_data;
    m_shape = other.m_shape;
    other.m_data = nullptr;
    other.m_shape = {0ul, 0ul};
  };

  T *m_data{};
  Shape m_shape{};
  Deleter m_deleter{};
  Acts::Device m_device{};
};

/// Element-wise sigmoid function for float cpu tensors
/// @param tensor The tensor to apply the sigmoid function to
/// @param execContext The execution context
void sigmoid(Tensor<float> &tensor, std::optional<cudaStream_t> stream = {});

/// Apply a score cut to the tensor and return a new tensor with the values that
/// satisfy the cut
/// @param tensor The tensor to apply the cut to
/// @param cut The cut value
/// @param execContext The execution context
std::pair<Tensor<float>, Tensor<std::int64_t>> applyScoreCut(
    const Tensor<float> &scores, const Tensor<std::int64_t> &edgeIndex,
    float cut, std::optional<cudaStream_t> stream = {});

}  // namespace Acts
