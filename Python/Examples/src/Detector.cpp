// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorCommons/Detector.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsExamples/DetectorCommons/StructureSelector.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/GenericDetector/AlignedGenericDetector.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;

namespace ActsPython {

void addDetector(py::module& mex) {
  {
    py::class_<Detector, std::shared_ptr<Detector>>(mex, "DetectorBase")
        .def("nominalGeometryContext", &Detector::nominalGeometryContext)
        .def("trackingGeometry", &Detector::trackingGeometry)
        .def("contextDecorators", &Detector::contextDecorators)
        .def("__enter__",
             [](const std::shared_ptr<Detector>& self) { return self; })
        .def("__exit__",
             [](std::shared_ptr<Detector>& self,
                const std::optional<py::object>&,
                const std::optional<py::object>&,
                const std::optional<py::object>&) { self.reset(); });
  }

  {
    py::class_<StructureSelector, std::shared_ptr<StructureSelector>>(
        mex, "StructureSelector")
        .def(py::init<std::shared_ptr<const TrackingGeometry>>())
        .def("selectSurfaces", &StructureSelector::selectSurfaces)
        .def("selectedTransforms", &StructureSelector::selectedTransforms);
  }

  {
    auto d =
        py::class_<GenericDetector, Detector, std::shared_ptr<GenericDetector>>(
            mex, "GenericDetector")
            .def(py::init<const GenericDetector::Config&>());

    auto c = py::class_<GenericDetector::Config>(d, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, buildLevel, logLevel, surfaceLogLevel, layerLogLevel,
                       volumeLogLevel, buildProto, materialDecorator, gen3,
                       graphvizFile);
  }

  {
    auto ad = py::class_<AlignedGenericDetector, GenericDetector,
                         std::shared_ptr<AlignedGenericDetector>>(
                  mex, "AlignedGenericDetector")
                  .def(py::init<const GenericDetector::Config&>());
  }

  {
    auto d =
        py::class_<TelescopeDetector, Detector,
                   std::shared_ptr<TelescopeDetector>>(mex, "TelescopeDetector")
            .def(py::init<const TelescopeDetector::Config&>());

    auto c =
        py::class_<TelescopeDetector::Config>(d, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, positions, stereos, offsets, bounds, thickness,
                       surfaceType, binValue, materialDecorator, logLevel);
  }
}

}  // namespace ActsPython
