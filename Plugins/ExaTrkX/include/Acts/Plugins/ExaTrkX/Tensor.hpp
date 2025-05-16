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
#include <optional>
#include <ostream>

#ifdef ACTS_EXATRKX_WITH_CUDA
typedef __device_builtin__ struct CUstream_st *cudaStream_t;
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

namespace detail {

/// This class implements the memory management for the Acts::Tensor
class TensorMemoryImpl {
 public:
  using Deleter = std::function<void(void *)>;

  TensorMemoryImpl(std::size_t nbytes, const ExecutionContext &execContext);

  TensorMemoryImpl(const TensorMemoryImpl &) = delete;
  TensorMemoryImpl(TensorMemoryImpl &&) noexcept;

  TensorMemoryImpl &operator=(const TensorMemoryImpl &) = delete;
  TensorMemoryImpl &operator=(TensorMemoryImpl &&) noexcept;

  ~TensorMemoryImpl();

  TensorMemoryImpl clone(std::size_t nbytes, const ExecutionContext &to) const;

  void *data() { return m_ptr; }
  const void *data() const { return m_ptr; }

  Acts::Device device() const { return m_device; }

 private:
  void moveConstruct(TensorMemoryImpl &&) noexcept;

  void *m_ptr{};
  Deleter m_deleter;
  Acts::Device m_device{};
};
}  // namespace detail

/// This is a very small, limited class that models a 2D tensor of arbitrary
/// type. It is move-only, and only possible to create via static factory
/// functions to ensure lifetime management.
/// This on purpose does not implement operations such as clone/to-host/to-cuda
template <typename T>
class Tensor {
 public:
  using Shape = std::array<std::size_t, 2>;

  static Tensor Create(Shape shape, const ExecutionContext &execContext) {
    detail::TensorMemoryImpl memory(shape[0] * shape[1] * sizeof(T),
                                    execContext);
    return Tensor(std::move(memory), shape);
  }

  /// Clone the tensor, copying the data to the new device
  /// @param to The {device, stream} to clone to
  /// @note This is a always a deep copy, even if the source and destination are the
  /// same device
  Tensor clone(const ExecutionContext &to) const {
    auto clonedMemory = m_memory.clone(nbytes(), to);
    return Tensor(std::move(clonedMemory), m_shape);
  }

  T *data() { return static_cast<T *>(m_memory.data()); }
  const T *data() const { return static_cast<const T *>(m_memory.data()); }
  Shape shape() const { return m_shape; }
  std::size_t size() const { return m_shape[0] * m_shape[1]; }
  std::size_t nbytes() const { return size() * sizeof(T); }
  Acts::Device device() const { return m_memory.device(); }

 private:
  Tensor(detail::TensorMemoryImpl memory, Shape shape)
      : m_shape(shape), m_memory(std::move(memory)) {}

  Shape m_shape{};
  detail::TensorMemoryImpl m_memory;
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
