// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/Cuda/CudaTracccBackend.hpp"

#include <traccc/cuda/utils/make_magnetic_field.hpp>
#include <traccc/cuda/utils/stream_wrapper.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/stream_wrapper.hpp>

namespace ActsExamples {

struct CudaTracccBackend::Impl {
  explicit Impl(int device)
      : deviceMr(device),
        mr(deviceMr, &hostMr),
        vecmemStream(device),
        stream(vecmemStream.stream()),
        copy(vecmemStream) {}

  vecmem::cuda::host_memory_resource hostMr;
  vecmem::cuda::device_memory_resource deviceMr;
  traccc::memory_resource mr;
  vecmem::cuda::stream_wrapper vecmemStream;
  traccc::cuda::stream_wrapper stream;
  mutable vecmem::cuda::async_copy copy;
};

CudaTracccBackend::CudaTracccBackend(const Config& cfg,
                                     std::unique_ptr<const Acts::Logger> _logger)
    : m_cfg(cfg),
      m_logger(_logger != nullptr
                   ? std::move(_logger)
                   : Acts::getDefaultLogger("CudaTracccBackend",
                                            Acts::Logging::INFO)),
      m_impl(std::make_unique<Impl>(cfg.device)) {
  ACTS_INFO("Set up traccc CUDA backend on device " << m_impl->stream.device());
}

CudaTracccBackend::~CudaTracccBackend() = default;

const traccc::memory_resource& CudaTracccBackend::memoryResource() const {
  return m_impl->mr;
}

vecmem::copy& CudaTracccBackend::copy() const {
  return m_impl->copy;
}

void CudaTracccBackend::synchronize() const {
  m_impl->stream.synchronize();
}

traccc::magnetic_field CudaTracccBackend::makeMagneticField(
    const traccc::magnetic_field& hostField) const {
  return traccc::cuda::make_magnetic_field(hostField);
}

}  // namespace ActsExamples
