// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Gnn/Stages.hpp"
#include "ActsPlugins/Gnn/Tensor.hpp"

#include <cstddef>
#include <vector>

namespace ActsPlugins::detail {

/// Every operation in this plugin that needs a GPU, declared without naming a
/// CUDA type: `cudaStream_t` here is the forward-declared `CUstream_st *` from
/// Tensor.hpp, which is the type the CUDA runtime uses too.
///
/// Exactly one implementation is linked -- the `.cu` translation units when the
/// plugin is built with CUDA, DeviceOpsNoCuda.cpp otherwise -- so callers below
/// dispatch on the runtime device and need no preprocessor branch. Both
/// implementations are built by CI, which is what keeps them in step.
///
/// Reaching any of these without CUDA means a CUDA tensor was constructed
/// first, which cudaCreateTensorMemory refuses, so the no-CUDA definitions are
/// unreachable rather than merely unsupported.

/// Allocate device memory for a tensor. @p ctx must carry a stream.
TensorPtr cudaCreateTensorMemory(std::size_t nbytes,
                                 const ExecutionContext &ctx);

/// Copy @p nbytes between host and device in either direction. @p to must carry
/// a stream.
void cudaCopyTensorMemory(void *dst, const void *src, std::size_t nbytes,
                          Device from, const ExecutionContext &to);

/// Apply the logistic function to @p tensor in place.
void cudaSigmoid(Tensor<float> &tensor, cudaStream_t stream);

/// Element-wise `scores > cut`.
Tensor<bool> cudaScoreMask(const Tensor<float> &scores, float cut,
                           cudaStream_t stream);

/// Gather the rows of @p tensor selected by @p mask.
template <typename T>
Tensor<T> cudaSelectRows(const Tensor<T> &tensor, const Tensor<bool> &mask,
                         const ExecutionContext &execContext);

/// Gather the columns of @p tensor selected by @p mask.
template <typename T>
Tensor<T> cudaSelectCols(const Tensor<T> &tensor, const Tensor<bool> &mask,
                         const ExecutionContext &execContext);

/// Scale each column of @p src by the matching entry of @p scales.
template <typename T>
Tensor<T> cudaMulPerColumn(const Tensor<T> &src, const Tensor<T> &scales,
                           const ExecutionContext &execContext);

/// Drop nodes no edge refers to, renumbering the edge index to match.
PipelineTensors cudaRemoveUnusedNodes(PipelineTensors &&tensors,
                                      std::vector<int> &spacePointIds,
                                      const ExecutionContext &execCtx);

}  // namespace ActsPlugins::detail
