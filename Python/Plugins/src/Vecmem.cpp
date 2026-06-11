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
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsVecmem, vecmem) {
  {
    py::class_<vecmem::memory_resource,
               std::shared_ptr<vecmem::memory_resource>>(vecmem,
                                                         "MemoryResource");

    py::class_<vecmem::host_memory_resource, vecmem::memory_resource,
               std::shared_ptr<vecmem::host_memory_resource>>(
        vecmem, "HostMemoryResource")
        .def(py::init<>());
  }
}
