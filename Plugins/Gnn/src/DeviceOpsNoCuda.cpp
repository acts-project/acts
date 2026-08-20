// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cstdint>

#include "CudaUnavailable.hpp"
#include "DeviceOps.hpp"

namespace ActsPlugins::detail {

TensorPtr cudaCreateTensorMemory(std::size_t /*nbytes*/,
                                 const ExecutionContext & /*ctx*/) {
  throwNoCudaSupport("create CUDA tensor");
}

void cudaCopyTensorMemory(void * /*dst*/, const void * /*src*/,
                          std::size_t /*nbytes*/, Device /*from*/,
                          const ExecutionContext & /*to*/) {
  throwNoCudaSupport("copy CUDA tensor memory");
}

void cudaSigmoid(Tensor<float> & /*tensor*/, cudaStream_t /*stream*/) {
  throwNoCudaSupport("apply sigmoid to CUDA tensor");
}

Tensor<bool> cudaScoreMask(const Tensor<float> & /*scores*/, float /*cut*/,
                           cudaStream_t /*stream*/) {
  throwNoCudaSupport("compute score mask on CUDA tensor");
}

template <typename T>
Tensor<T> cudaSelectRows(const Tensor<T> & /*tensor*/,
                         const Tensor<bool> & /*mask*/,
                         const ExecutionContext & /*execContext*/) {
  throwNoCudaSupport("selectRows on CUDA tensor");
}

template <typename T>
Tensor<T> cudaSelectCols(const Tensor<T> & /*tensor*/,
                         const Tensor<bool> & /*mask*/,
                         const ExecutionContext & /*execContext*/) {
  throwNoCudaSupport("selectCols on CUDA tensor");
}

template <typename T>
Tensor<T> cudaMulPerColumn(const Tensor<T> & /*src*/,
                           const Tensor<T> & /*scales*/,
                           const ExecutionContext & /*execContext*/) {
  throwNoCudaSupport("apply feature scales on CUDA tensor");
}

PipelineTensors cudaRemoveUnusedNodes(PipelineTensors && /*tensors*/,
                                      std::vector<int> & /*spacePointIds*/,
                                      const ExecutionContext & /*execCtx*/) {
  throwNoCudaSupport("removeUnusedNodes on CUDA tensor");
}

// The same set Tensor.cu instantiates. A mismatch here is a link error in the
// CUDA-off build, which is why that build is part of CI.
template Tensor<float> cudaSelectRows(const Tensor<float> &,
                                      const Tensor<bool> &,
                                      const ExecutionContext &);
template Tensor<std::int64_t> cudaSelectRows(const Tensor<std::int64_t> &,
                                             const Tensor<bool> &,
                                             const ExecutionContext &);
template Tensor<float> cudaSelectCols(const Tensor<float> &,
                                      const Tensor<bool> &,
                                      const ExecutionContext &);
template Tensor<std::int64_t> cudaSelectCols(const Tensor<std::int64_t> &,
                                             const Tensor<bool> &,
                                             const ExecutionContext &);
template Tensor<float> cudaMulPerColumn(const Tensor<float> &,
                                        const Tensor<float> &,
                                        const ExecutionContext &);
template Tensor<std::int64_t> cudaMulPerColumn(const Tensor<std::int64_t> &,
                                               const Tensor<std::int64_t> &,
                                               const ExecutionContext &);

}  // namespace ActsPlugins::detail
