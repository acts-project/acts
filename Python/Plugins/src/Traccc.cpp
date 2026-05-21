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

#if defined(ACTS_ENABLE_CUDA)
#include <traccc/cuda/utils/stream.hpp>

void bind_traccc_cuda(py::module_& m){
    py::class_<traccc::cuda::stream,
               std::shared_ptr<traccc::cuda::stream>>(m, "CudaStream")
        .def(py::init<>());

}
#endif

namespace py = pybind11;
using namespace pybind11::literals;

void bind_traccc_host(py::module_& m){
    // ---- traccc types ----
    py::class_<traccc::memory_resource,
               std::shared_ptr<traccc::memory_resource>>(m,
                                                         "TracccMemoryResource")
        .def(py::init([](vecmem::memory_resource& main,
                         vecmem::memory_resource* host) {
               return std::make_shared<traccc::memory_resource>(main, host);
             }),
             py::arg("main"), py::arg("host") = nullptr,
             py::keep_alive<1, 2>(),   // keep main alive
             py::keep_alive<1, 3>());  // keep host alive
}

PYBIND11_MODULE(ActsPluginsPythonBindingsTraccc, traccc) {

    bind_traccc_host(traccc);
#if defined(ACTS_ENABLE_CUDA)
    bind_traccc_cuda(traccc);
#endif
    
}