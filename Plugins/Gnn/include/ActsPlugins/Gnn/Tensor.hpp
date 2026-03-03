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

namespace ActsPlugins {

/// @addtogroup gnn_plugin
/// @{

/// A simple device description struct
struct Device {
  /// Device type enumeration
  enum class Type { eCPU, eCUDA };
  /// Device type (CPU or CUDA)
  Type type = Type::eCPU;
  /// Device index for multi-GPU systems
  std::size_t index = 0;

  /// @brief Create a CPU device descriptor
  /// @return Device object configured for CPU execution
  static Device Cpu() { return {Type::eCPU, 0}; }
  /// @brief Create a CUDA device descriptor
  /// @param index GPU device index for multi-GPU systems (default: 0)
  /// @return Device object configured for CUDA execution on specified GPU
  static Device Cuda(std::size_t index = 0) { return {Type::eCUDA, index}; }

  /// @brief Check if device is configured for CPU execution
  /// @return True if device type is CPU, false otherwise
  bool isCpu() const { return type == Type::eCPU; }
  /// @brief Check if device is configured for CUDA execution
  /// @return True if device type is CUDA, false otherwise
  bool isCuda() const { return type == Type::eCUDA; }

  /// @brief Compare two device descriptors for equality
  /// @return True if both devices have same type and index
  bool operator==(const Device &) const = default;
  /// @brief Compare two device descriptors for equality
  /// @param other Device to compare against
  /// @return True if both devices have same type and index
};

/// Stream operator for Device
/// @param os Output stream
/// @param device Device to output
/// @return Reference to output stream
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
  /// Target device for execution
  Device device{Device::Type::eCPU};
  /// CUDA stream for asynchronous execution
  std::optional<cudaStream_t> stream;
};

/// @cond
namespace detail {

using TensorDeleter = std::function<void(void *)>;
using TensorPtr = std::unique_ptr<void, TensorDeleter>;

TensorPtr createTensorMemory(std::size_t nbytes, const ExecutionContext &ctx);
TensorPtr cloneTensorMemory(const TensorPtr &ptrFrom, std::size_t nbytes,
                            Device devFrom, const ExecutionContext &ctxTo);

}  // namespace detail
/// @endcond

/// This is a very small, limited class that models a 2D tensor of arbitrary
/// type. It is move-only, and only possible to create via static factory
/// functions to ensure lifetime management.
/// This on purpose does not implement operations such as clone/to-host/to-cuda
template <Acts::Concepts::arithmetic T>
class Tensor {
 public:
  /// Type alias for tensor shape as 2D dimensions array
  using Shape = std::array<std::size_t, 2>;

  /// @brief Create a new tensor with specified shape and execution context
  /// @param shape 2D tensor dimensions [rows, columns]
  /// @param execContext Execution context specifying device and optional CUDA stream
  /// @return Newly created tensor with allocated memory on specified device
  static Tensor Create(Shape shape, const ExecutionContext &execContext) {
    auto ptr = detail::createTensorMemory(shape[0] * shape[1] * sizeof(T),
                                          execContext);
    return Tensor(shape, std::move(ptr), execContext);
  }

  /// Clone the tensor, copying the data to the new device
  /// @param to The {device, stream} to clone to
  /// @note This is a always a deep copy, even if the source and destination are the
  /// same device
  /// @return New tensor on target device with copied data
  Tensor clone(const ExecutionContext &to) const {
    auto clonedPtr = detail::cloneTensorMemory(m_ptr, nbytes(), m_device, to);
    return Tensor(m_shape, std::move(clonedPtr), to);
  }

  /// Get the non-const data pointer
  /// @return Mutable pointer to tensor data
  T *data() { return static_cast<T *>(m_ptr.get()); }

  /// Get the const data pointer
  /// @return Const pointer to tensor data
  const T *data() const { return static_cast<const T *>(m_ptr.get()); }

  /// Get the shape of the tensor
  /// @return Array containing tensor dimensions [rows, columns]
  Shape shape() const { return m_shape; }

  /// Get the number of elements in the tensor
  /// @return Total number of elements in the tensor
  std::size_t size() const { return m_shape[0] * m_shape[1]; }

  /// Get the number of bytes of the tensor
  /// @return Total memory size of tensor data in bytes
  std::size_t nbytes() const { return size() * sizeof(T); }

  /// Get the device of the tensor
  /// @return Device where tensor data is stored
  Device device() const { return m_device; }

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
/// @return Pair of filtered score and edge index tensors containing only edges above threshold
std::pair<Tensor<float>, Tensor<std::int64_t>> applyScoreCut(
    const Tensor<float> &scores, const Tensor<std::int64_t> &edgeIndex,
    float cut, std::optional<cudaStream_t> stream = {});

/// Apply a limit on the number of edges consistently on edgeIndex and
/// edgeFeatures.
/// @param edgeIndex The edge index tensor
/// @param edgeFeatures The edge feature tensor
/// @param maxEdges The edge limit to apply
/// @param stream The stream to use for operation in case of CUDA
/// @return Pair of limited edge index tensor and optional limited edge features tensor
std::pair<Tensor<std::int64_t>, std::optional<Tensor<float>>> applyEdgeLimit(
    const Tensor<std::int64_t> &edgeIndex,
    const std::optional<Tensor<float>> &edgeFeatures, std::size_t maxEdges,
    std::optional<cudaStream_t> stream);

/// @}
}  // namespace ActsPlugins
