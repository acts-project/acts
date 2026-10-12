// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Nodd/DD4hepNoddDetector.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsNodd, m) {
  // Register the inherited Detector and DD4hep configuration types first.
  py::module_::import("acts.examples.dd4hep");

  auto detector =
      py::class_<DD4hepNoddDetector, DD4hepDetectorBase,
                 std::shared_ptr<DD4hepNoddDetector>>(m, "DD4hepNoddDetector")
          .def(py::init<const DD4hepNoddDetector::Config&>())
          .def("geometryReport", &DD4hepNoddDetector::geometryReport)
          .def_property_readonly("config", &DD4hepNoddDetector::config,
                                 py::return_value_policy::reference_internal);

  auto config =
      py::class_<DD4hepNoddDetector::Config, DD4hepDetectorBase::Config>(
          detector, "Config")
          .def(py::init<>());
  ACTS_PYTHON_STRUCT(config, gen3, navigationEnvelopeMm);
  patchKwargsConstructor(config);
}
