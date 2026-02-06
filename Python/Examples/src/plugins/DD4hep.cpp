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
    auto base =
        py::class_<DD4hepDetectorBase, Detector,
                   std::shared_ptr<DD4hepDetectorBase>>(m, "DD4hepDetectorBase")
            .def_property_readonly("field", &DD4hepDetectorBase::field);
    auto c = py::class_<DD4hepDetectorBase::Config>(base, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, logLevel, dd4hepLogLevel, xmlFileNames, name);
    patchKwargsConstructor(c);
  }

  {
    auto f = py::class_<DD4hepDetector, DD4hepDetectorBase,
                        std::shared_ptr<DD4hepDetector>>(m, "DD4hepDetector")
                 .def(py::init<const DD4hepDetector::Config&>());

    auto c = py::class_<DD4hepDetector::Config, DD4hepDetectorBase::Config>(
                 f, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, bTypePhi, bTypeR, bTypeZ, envelopeR, envelopeZ,
                       defaultLayerThickness, materialDecorator,
                       geometryIdentifierHook, detectorElementFactory);
    patchKwargsConstructor(c);

    m.def("alignedDD4hepDetectorElementFactory",
          &alignedDD4hepDetectorElementFactory);
  }

  {
    auto odd =
        py::class_<OpenDataDetector, DD4hepDetectorBase,
                   std::shared_ptr<OpenDataDetector>>(m, "OpenDataDetector")
            .def(py::init<const OpenDataDetector::Config&,
                          const Acts::GeometryContext&>(),
                 "config"_a, "gctx"_a);

    auto c = py::class_<OpenDataDetector::Config, DD4hepDetectorBase::Config>(
                 odd, "Config")
                 .def(py::init<>());
    // ACTS_PYTHON_STRUCT(c, );

    patchKwargsConstructor(c);
  }
}
