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

#include <traccc/cuda/utils/stream.hpp>
#include <traccc/utils/memory_resource.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsVecmem, vecmem) {
  auto cuda = vecmem.def_submodule("cuda", "CUDA backend bindings");

  py::class_<vecmem::cuda::host_memory_resource, vecmem::memory_resource,
             std::shared_ptr<vecmem::cuda::host_memory_resource>>(
      cuda, "host_memory_resource")
      .def(py::init<>());

  py::class_<vecmem::cuda::device_memory_resource, vecmem::memory_resource,
             std::shared_ptr<vecmem::cuda::device_memory_resource>>(
      cuda, "device_memory_resource")
      .def(py::init<>());

  py::class_<traccc::cuda::stream, std::shared_ptr<traccc::cuda::stream>>(
      cuda, "stream")
      .def(py::init<>());

  py::class_<vecmem::cuda::async_copy, vecmem::copy,
             std::shared_ptr<vecmem::cuda::async_copy>>(cuda, "async_copy")
      .def(py::init([](traccc::cuda::stream& s) {
             return std::make_shared<vecmem::cuda::async_copy>(s.cudaStream());
           }),
           py::arg("stream"), py::keep_alive<1, 2>());
}
