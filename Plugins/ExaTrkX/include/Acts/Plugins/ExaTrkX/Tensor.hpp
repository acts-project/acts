// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Concepts.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <memory>
#include <optional>
#include <ostream>

/// Forward declare cuda stream, to be able to use the header without cuda
struct CUstream_st;
using cudaStream_t = CUstream_st *;

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

using TensorDeleter = std::function<void(void *)>;
using TensorPtr = std::unique_ptr<void, TensorDeleter>;

TensorPtr createTensorMemory(std::size_t nbytes, const ExecutionContext &ctx);
TensorPtr cloneTensorMemory(const TensorPtr &ptrFrom, std::size_t nbytes,
                            Acts::Device devFrom,
                            const ExecutionContext &ctxTo);

}  // namespace detail

/// This is a very small, limited class that models a 2D tensor of arbitrary
/// type. It is move-only, and only possible to create via static factory
/// functions to ensure lifetime management.
/// This on purpose does not implement operations such as clone/to-host/to-cuda
template <Acts::Concepts::arithmetic T>
class Tensor {
 public:
  using Shape = std::array<std::size_t, 2>;

  static Tensor Create(Shape shape, const ExecutionContext &execContext) {
    auto ptr = detail::createTensorMemory(shape[0] * shape[1] * sizeof(T),
                                          execContext);
    return Tensor(shape, std::move(ptr), execContext);
  }

  /// Clone the tensor, copying the data to the new device
  /// @param to The {device, stream} to clone to
  /// @note This is a always a deep copy, even if the source and destination are the
  /// same device
  Tensor clone(const ExecutionContext &to) const {
    auto clonedPtr = detail::cloneTensorMemory(m_ptr, nbytes(), m_device, to);
    return Tensor(m_shape, std::move(clonedPtr), to);
  }

  /// Get the non-const data pointer
  T *data() { return static_cast<T *>(m_ptr.get()); }

  /// Get the const data pointer
  const T *data() const { return static_cast<const T *>(m_ptr.get()); }

  /// Get the shape of the tensor
  Shape shape() const { return m_shape; }

  /// Get the number of elements in the tensor
  std::size_t size() const { return m_shape[0] * m_shape[1]; }

  /// Get the number of bytes of the tensor
  std::size_t nbytes() const { return size() * sizeof(T); }

  /// Get the device of the tensor
  Acts::Device device() const { return m_device; }

 private:
  Tensor(Shape shape, detail::TensorPtr ptr, const ExecutionContext &ctx)
      : m_shape(shape), m_ptr(std::move(ptr)), m_device(ctx.device) {}

  Shape m_shape{};
  detail::TensorPtr m_ptr;
  Device m_device{};
};

/// Element-wise sigmoid function for float cpu tensors
/// @param tensor The tensor to apply the sigmoid function to
/// @param stream The stream to use for the operation in case of CUDA
void sigmoid(Tensor<float> &tensor, std::optional<cudaStream_t> stream = {});

/// Apply a score cut to the tensor and return a new tensor with the values that
/// satisfy the cut
/// @param scores The edge score tensor
/// @param edgeIndex The edge index tensor
/// @param cut The score cut value which edges to accept
/// @param stream The stream to use for the operation in case of CUDA
std::pair<Tensor<float>, Tensor<std::int64_t>> applyScoreCut(
    const Tensor<float> &scores, const Tensor<std::int64_t> &edgeIndex,
    float cut, std::optional<cudaStream_t> stream = {});

}  // namespace Acts
