// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// CUDA include(s).
#include <cuda_runtime_api.h>

// System include(s).
#include <memory>
#include <string>

namespace traccc::cuda::details {

/// Hough seeding peak filter network, run with ONNX Runtime
///
/// The model takes images [n, 32, 32, 1] and positions [n, 1] and writes
/// scores [n, 2]. All pointers are device memory aligned to 256 bytes.
///
class peak_filter {
 public:
  static constexpr unsigned int image_size = 32u * 32u;

  peak_filter(const std::string& model_path, cudaStream_t stream);
  ~peak_filter();

  void operator()(const float* images, const float* positions, unsigned int n,
                  float* scores) const;

 private:
  struct impl;
  std::unique_ptr<impl> m_impl;
};

}  // namespace traccc::cuda::details
