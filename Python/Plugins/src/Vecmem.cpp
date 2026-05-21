// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <string>

#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>

#if defined(ACTS_ENABLE_CUDA)

#include <traccc/cuda/utils/stream.hpp>
#include <traccc/utils/memory_resource.hpp>
#include <vecmem/copy/cuda/async_copy.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>

void bind_vecmem_cuda(py::module_& m) {
  py::class_<vecmem::cuda::host_memory_resource, vecmem::memory_resource,
             std::shared_ptr<vecmem::cuda::host_memory_resource>>(
      m, "CudaHostMemoryResource")
      .def(py::init<>());

  py::class_<vecmem::cuda::device_memory_resource, vecmem::memory_resource,
             std::shared_ptr<vecmem::cuda::device_memory_resource>>(
      m, "CudaDeviceMemoryResource")
      .def(py::init<>());

  py::class_<traccc::cuda::stream, std::shared_ptr<traccc::cuda::stream>>(
      m, "CudaStream")
      .def(py::init<>());

  py::class_<vecmem::cuda::async_copy, vecmem::copy,
             std::shared_ptr<vecmem::cuda::async_copy>>(m, "CudaAsyncCopy")
      .def(py::init([](traccc::cuda::stream& s) {
             return std::make_shared<vecmem::cuda::async_copy>(s.cudaStream());
           }),
           py::arg("stream"), py::keep_alive<1, 2>());

  py::class_<traccc::memory_resource, std::shared_ptr<traccc::memory_resource>>(
      m, "TracccMemoryResource")
      .def(py::init([](vecmem::memory_resource& main,
                       vecmem::memory_resource* host) {
             return std::make_shared<traccc::memory_resource>(main, host);
           }),
           py::arg("main"), py::arg("host") = nullptr, py::keep_alive<1, 2>(),
           py::keep_alive<1, 3>());
}
#endif

namespace py = pybind11;
using namespace pybind11::literals;

void bind_vecmem_host(py::module_& m) {
  // HostMemory composite to be passed to algorithms expecting both copy and
  // memory resource
  py::class_<vecmem::memory_resource, std::shared_ptr<vecmem::memory_resource>>(
      m, "MemoryResource");
  py::class_<vecmem::copy, std::shared_ptr<vecmem::copy>>(m, "Copy");
  py::class_<vecmem::host_memory_resource, vecmem::memory_resource,
             std::shared_ptr<vecmem::host_memory_resource>>(
      m, "HostMemoryResource")
      .def(py::init<>());
  py::class_<vecmem::copy, std::shared_ptr<vecmem::copy>>(m, "HostCopy")
      .def(py::init<>());
}

PYBIND11_MODULE(ActsPluginsPythonBindingsVecmem, m) {
  bind_vecmem_host(m);

// registers device memoory composites corresponding to the backends enabled
// during build
#if defined(ACTS_ENABLE_CUDA)
  bind_vecmem_cuda(m);
#endif
}
