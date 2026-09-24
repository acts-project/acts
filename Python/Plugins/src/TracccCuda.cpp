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

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsTraccc, traccc) {
  auto cuda = traccc.def_submodule("cuda", "CUDA backend bindings");

  py::class_<traccc::cuda::stream, std::shared_ptr<traccc::cuda::stream>>(
      cuda, "stream")
      .def(py::init<>());
}
