// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DD4hepDetector/AlignedDD4hepDetectorElement.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/DD4hepDetector/OpenDataDetector.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <utility>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Acts;
using namespace ActsPlugins;
using namespace pybind11::literals;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsDD4hep, m) {
  {
    py::class_<AlignedDD4hepDetectorElement, DD4hepDetectorElement,
               std::shared_ptr<AlignedDD4hepDetectorElement>>(
        m, "AlignedDD4hepDetectorElement");
  }

  {
    auto f =
        py::class_<DD4hepDetector, Detector, std::shared_ptr<DD4hepDetector>>(
            m, "DD4hepDetector")
            .def(py::init<const DD4hepDetector::Config&>())
            .def_property_readonly("field", &DD4hepDetector::field);

    auto c = py::class_<DD4hepDetector::Config>(f, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, logLevel, dd4hepLogLevel, xmlFileNames, name,
                       bTypePhi, bTypeR, bTypeZ, envelopeR, envelopeZ,
                       defaultLayerThickness, materialDecorator,
                       geometryIdentifierHook, detectorElementFactory);
    patchKwargsConstructor(c);

    m.def("alignedDD4hepDetectorElementFactory",
          &alignedDD4hepDetectorElementFactory);
  }
}
