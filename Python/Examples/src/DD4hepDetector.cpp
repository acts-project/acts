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
#include "ActsPython/Utilities/Macros.hpp"
#include "ActsPython/Utilities/Patchers.hpp"

#include <DD4hep/Fields.h>

#include <memory>
#include <utility>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace ActsExamples;

namespace ActsPython {

/// This adds the TGeoDetector to the examples module
/// @param mex the examples module
void addDD4hepDetector(py::module_& mex) {
  auto dd4hep = mex.def_submodule("dd4hep");

  {

    py::class_<ActsExamples::AlignedDD4hepDetectorElement,
               Acts::DD4hepDetectorElement,
               std::shared_ptr<ActsExamples::AlignedDD4hepDetectorElement>>(
        dd4hep, "AlignedDD4hepDetectorElement");
  }

  {
    auto f =
        py::class_<DD4hepDetector, Detector, std::shared_ptr<DD4hepDetector>>(
            dd4hep, "DD4hepDetector")
            .def(py::init<const DD4hepDetector::Config&>())
            .def_property_readonly("field", &DD4hepDetector::field);

    auto c = py::class_<DD4hepDetector::Config>(f, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(
        c, logLevel, dd4hepLogLevel, xmlFileNames, name, bTypePhi, bTypeR,
        bTypeZ, envelopeR, envelopeZ, defaultLayerThickness, materialDecorator,
        alignmentDecorator, geometryIdentifierHook, detectorElementFactory);
    patchKwargsConstructor(c);

    dd4hep.def("alignedDD4hepDetectorElementFactory",
          &ActsExamples::alignedDD4hepDetectorElementFactory);
  }

  {
    py::class_<Acts::DD4hepFieldAdapter, Acts::MagneticFieldProvider,
               std::shared_ptr<Acts::DD4hepFieldAdapter>>(dd4hep,
                                                          "DD4hepFieldAdapter");
  }

}

} // namespace ActsPython
