// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/DetectorCommons/AlignmentDecorator.hpp"
#include "ActsExamples/DetectorCommons/AlignmentGenerator.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;

namespace Acts::Python {

void addAlignment(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    auto ad = py::class_<AlignmentDecorator, IContextDecorator,
                         std::shared_ptr<AlignmentDecorator>>(
                  mex, "AlignmentDecorator")
                  .def(py::init<const AlignmentDecorator::Config&,
                                Acts::Logging::Level>())
                  .def("decorate", &AlignmentDecorator::decorate)
                  .def("name", &AlignmentDecorator::name);

    auto c =
        py::class_<AlignmentDecorator::Config>(ad, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, iovStores, nominalStore, garbageCollection,
                       gcInterval, iovGenerators);
  }

  {
    py::class_<IAlignmentStore, std::shared_ptr<IAlignmentStore>>(
        mex, "IAlignmentStore");
  }

  {
    py::class_<GeoIdAlignmentStore, IAlignmentStore,
               std::shared_ptr<GeoIdAlignmentStore>>(mex, "GeoIdAlignmentStore")
        .def(py::init<const std::unordered_map<Acts::GeometryIdentifier,
                                               Acts::Transform3>&>());
  }

  {
    py::class_<AlignmentGenerator::Nominal>(mex, "AlignmentGeneratorNominal")
        .def(py::init<>())
        .def("__call__", &AlignmentGenerator::Nominal::operator());
  }

  {
    py::class_<AlignmentGenerator::GlobalShift>(mex,
                                                "AlignmentGeneratorGlobalShift")
        .def(py::init<>())
        .def_readwrite("shift", &AlignmentGenerator::GlobalShift::shift)
        .def_readwrite("randomize", &AlignmentGenerator::GlobalShift::randomize)
        .def("__call__", &AlignmentGenerator::GlobalShift::operator());
  }

  {
    py::class_<AlignmentGenerator::GlobalRotation>(
        mex, "AlignmentGeneratorGlobalRotation")
        .def(py::init<>())
        .def_readwrite("axis", &AlignmentGenerator::GlobalRotation::axis)
        .def_readwrite("angle", &AlignmentGenerator::GlobalRotation::angle)
        .def_readwrite("randomize",
                       &AlignmentGenerator::GlobalRotation::randomize)
        .def("__call__", &AlignmentGenerator::GlobalRotation::operator());
  }

  {
    py::class_<AlignmentGenerator::LocalRotation>(
        mex, "AlignmentGeneratorLocalRotation")
        .def(py::init<>())
        .def_readwrite("axis", &AlignmentGenerator::LocalRotation::axis)
        .def_readwrite("angle", &AlignmentGenerator::LocalRotation::angle)
        .def_readwrite("randomize",
                       &AlignmentGenerator::LocalRotation::randomize)
        .def("__call__", &AlignmentGenerator::LocalRotation::operator());
  }

  {
    py::class_<AlignmentGenerator::LocalShift>(mex,
                                               "AlignmentGeneratorLocalShift")
        .def(py::init<>())
        .def_readwrite("axisDirection",
                       &AlignmentGenerator::LocalShift::axisDirection)
        .def_readwrite("shift", &AlignmentGenerator::LocalShift::shift)
        .def_readwrite("randomize", &AlignmentGenerator::LocalShift::randomize)
        .def("__call__", &AlignmentGenerator::LocalShift::operator());
  }
}

}  // namespace Acts::Python
