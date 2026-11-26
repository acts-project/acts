// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "peak_filter.hpp"

#include "../utils/cuda_error_handling.hpp"

// ONNX Runtime include(s).
#include <onnxruntime_cxx_api.h>

// System include(s).
#include <algorithm>
#include <array>
#include <cstdint>

namespace traccc::cuda::details {
namespace {

/// Largest batch per network call, bounds the engine memory
constexpr unsigned int max_batch = 1024u;

std::string profile_shape(const unsigned int batch) {
  const std::string n = std::to_string(batch);
  return "image:" + n + "x32x32x1,positions:" + n + "x1";
}

}  // namespace

struct peak_filter::impl {
  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "traccc hough peak filter"};
  Ort::Session session{nullptr};
  Ort::MemoryInfo memory{nullptr};
};

peak_filter::peak_filter(const std::string& model_path, cudaStream_t stream)
    : m_impl(std::make_unique<impl>()) {
  int device = 0;
  TRACCC_CUDA_ERROR_CHECK(cudaGetDevice(&device));

  const std::string device_id = std::to_string(device);
  const std::string min_shape = profile_shape(1u);
  const std::string max_shape = profile_shape(max_batch);
  const std::array keys{"device_id", "trt_int8_enable",
                        "trt_profile_min_shapes", "trt_profile_opt_shapes",
                        "trt_profile_max_shapes"};
  const std::array values{device_id.c_str(), "1", min_shape.c_str(),
                          max_shape.c_str(), max_shape.c_str()};

  const auto& api = Ort::GetApi();
  OrtTensorRTProviderOptionsV2* trt_raw = nullptr;
  Ort::ThrowOnError(api.CreateTensorRTProviderOptions(&trt_raw));
  const std::unique_ptr<OrtTensorRTProviderOptionsV2,
                        decltype(api.ReleaseTensorRTProviderOptions)>
      trt_options(trt_raw, api.ReleaseTensorRTProviderOptions);
  Ort::ThrowOnError(api.UpdateTensorRTProviderOptions(
      trt_options.get(), keys.data(), values.data(), keys.size()));
  Ort::ThrowOnError(api.UpdateTensorRTProviderOptionsWithValue(
      trt_options.get(), "user_compute_stream", stream));

  Ort::SessionOptions options;
  options.AppendExecutionProvider_TensorRT_V2(*trt_options);

  m_impl->session = Ort::Session(m_impl->env, model_path.c_str(), options);
  m_impl->memory =
      Ort::MemoryInfo("Cuda", OrtDeviceAllocator, device, OrtMemTypeDefault);
}

peak_filter::~peak_filter() = default;

void peak_filter::operator()(const float* images, const float* positions,
                             const unsigned int n, float* scores) const {
  for (unsigned int first = 0u; first < n; first += max_batch) {
    const unsigned int count = std::min(max_batch, n - first);
    const std::array<std::int64_t, 4> image_shape{count, 32, 32, 1};
    const std::array<std::int64_t, 2> position_shape{count, 1};
    const std::array<std::int64_t, 2> score_shape{count, 2};

    Ort::Value image = Ort::Value::CreateTensor<float>(
        m_impl->memory, const_cast<float*>(images + first * image_size),
        count * image_size, image_shape.data(), image_shape.size());
    Ort::Value position = Ort::Value::CreateTensor<float>(
        m_impl->memory, const_cast<float*>(positions + first), count,
        position_shape.data(), position_shape.size());
    Ort::Value score = Ort::Value::CreateTensor<float>(
        m_impl->memory, scores + 2u * first, 2u * count, score_shape.data(),
        score_shape.size());

    Ort::IoBinding binding(m_impl->session);
    binding.BindInput("image", image);
    binding.BindInput("positions", position);
    binding.BindOutput("out", score);
    m_impl->session.Run(Ort::RunOptions{nullptr}, binding);
  }
}

}  // namespace traccc::cuda::details
