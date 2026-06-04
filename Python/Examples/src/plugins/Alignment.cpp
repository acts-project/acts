// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Alignment.hpp"

#include "ActsAlignment/Geometry/AlignableStructure.hpp"
#include "ActsAlignment/Geometry/AlignmentHierarchy.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"
#include "ActsExamples/DetectorCommons/AlignmentDecorator.hpp"
#include "ActsExamples/DetectorCommons/AlignmentGenerator.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsAlignment;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsAlignment, m) {
  {
    py::enum_<AlignmentMask>(m, "AlignmentMask")
        .value("None_", AlignmentMask::None)
        .value("Center0", AlignmentMask::Center0)
        .value("Center1", AlignmentMask::Center1)
        .value("Center2", AlignmentMask::Center2)
        .value("Rotation0", AlignmentMask::Rotation0)
        .value("Rotation1", AlignmentMask::Rotation1)
        .value("Rotation2", AlignmentMask::Rotation2)
        .value("All", AlignmentMask::All)
        .def(py::self | py::self)
        .def(py::self & py::self)
        .def(py::self ^ py::self)
        .def(~py::self);
  }

  {
    py::enum_<Acts::AlignmentIndices>(m, "AlignmentIndices")
        .value("eAlignmentCenter0", Acts::eAlignmentCenter0)
        .value("eAlignmentCenter1", Acts::eAlignmentCenter1)
        .value("eAlignmentCenter2", Acts::eAlignmentCenter2)
        .value("eAlignmentRotation0", Acts::eAlignmentRotation0)
        .value("eAlignmentRotation1", Acts::eAlignmentRotation1)
        .value("eAlignmentRotation2", Acts::eAlignmentRotation2)
        .value("eAlignmentSize", Acts::eAlignmentSize)
        .export_values();
  }

  {
    auto ad =
        py::class_<AlignmentDecorator, IContextDecorator,
                   std::shared_ptr<AlignmentDecorator>>(m, "AlignmentDecorator")
            .def(py::init<const AlignmentDecorator::Config&, Logging::Level>())
            .def("decorate", &AlignmentDecorator::decorate)
            .def("name", &AlignmentDecorator::name);

    auto c =
        py::class_<AlignmentDecorator::Config>(ad, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, iovStores, nominalStore, garbageCollection,
                       gcInterval, iovGenerators);
  }

  {
    py::class_<IAlignmentStore, std::shared_ptr<IAlignmentStore>>(
        m, "IAlignmentStore");
  }

  {
    py::class_<GeoIdAlignmentStore, IAlignmentStore,
               std::shared_ptr<GeoIdAlignmentStore>>(m, "GeoIdAlignmentStore")
        .def(py::init<
             const std::unordered_map<GeometryIdentifier, Transform3>&>());
  }

  {
    py::class_<AlignableStructure, std::shared_ptr<AlignableStructure>>(
        m, "AlignableStructure")
        .def(py::init<Acts::GeometryIdentifier>(), py::arg("id"))
        .def("addSurface", &AlignableStructure::addSurface, py::arg("surface"))
        .def("addChild", &AlignableStructure::addChild, py::arg("child"))
        .def_property_readonly("geometryId", &AlignableStructure::geometryId)
        .def_property_readonly(
            "surfaces",
            [](const AlignableStructure& s) {
              std::vector<const Acts::Surface*> result;
              for (const auto& srf : s.surfaces()) {
                result.push_back(&srf);
              }
              return result;
            },
            py::return_value_policy::reference_internal)
        .def_property_readonly(
            "children",
            static_cast<
                const std::vector<std::shared_ptr<AlignableStructure>>& (
                    AlignableStructure::*)() const>(
                &AlignableStructure::children),
            py::return_value_policy::reference_internal)
        .def_property(
            "alignmentMask",
            [](const AlignableStructure& s) { return s.alignmentMask(); },
            [](AlignableStructure& s, AlignmentMask m) {
              s.alignmentMask() = m;
            })
        .def("constraints",
             static_cast<
                 const std::unordered_map<Acts::AlignmentIndices, double>& (
                     AlignableStructure::*)() const>(
                 &AlignableStructure::constraints),
             py::return_value_policy::copy)
        .def(
            "setConstraint",
            [](AlignableStructure& s, Acts::AlignmentIndices idx,
               double variance) { s.constraints()[idx] = variance; },
            py::arg("index"), py::arg("variance"))
        .def(
            "clearConstraint",
            [](AlignableStructure& s, Acts::AlignmentIndices idx) {
              s.constraints().erase(idx);
            },
            py::arg("index"));
  }

  {
    auto hierarchy =
        py::class_<AlignmentHierarchy>(m, "AlignmentHierarchy")
            .def(py::init<
                     const std::vector<std::shared_ptr<AlignableStructure>>&>(),
                 py::arg("structures"))
            .def("validate", &AlignmentHierarchy::validate)
            .def("detectMaskConflicts",
                 &AlignmentHierarchy::detectMaskConflicts,
                 py::arg("moduleMask"), py::arg("floatingModules"))
            .def("structureFor",
                 static_cast<const AlignableStructure* (
                     AlignmentHierarchy::*)(const Acts::SurfacePlacementBase&)
                                 const>(&AlignmentHierarchy::structureFor),
                 py::arg("detElement"),
                 py::return_value_policy::reference_internal)
            .def("structureFor",
                 static_cast<const AlignableStructure* (
                     AlignmentHierarchy::*)(const Acts::Surface&) const>(
                     &AlignmentHierarchy::structureFor),
                 py::arg("surface"),
                 py::return_value_policy::reference_internal)
            .def_property_readonly("structures",
                                   &AlignmentHierarchy::structures,
                                   py::return_value_policy::reference_internal);

    py::class_<AlignmentHierarchy::ValidationResult>(hierarchy,
                                                     "ValidationResult")
        .def_readonly("overlapping",
                      &AlignmentHierarchy::ValidationResult::overlapping)
        .def("ok", &AlignmentHierarchy::ValidationResult::ok);

    py::class_<AlignmentHierarchy::MaskConflict>(hierarchy, "MaskConflict")
        .def_readonly("detElement",
                      &AlignmentHierarchy::MaskConflict::detElement)
        .def_readonly("structure", &AlignmentHierarchy::MaskConflict::structure)
        .def_readonly("conflictingBits",
                      &AlignmentHierarchy::MaskConflict::conflictingBits);
  }

  {
    py::class_<AlignmentGenerator::Nominal>(m, "AlignmentGeneratorNominal")
        .def(py::init<>())
        .def("__call__", &AlignmentGenerator::Nominal::operator());
  }

  {
    py::class_<AlignmentGenerator::GlobalShift>(m,
                                                "AlignmentGeneratorGlobalShift")
        .def(py::init<>())
        .def_readwrite("shift", &AlignmentGenerator::GlobalShift::shift)
        .def_readwrite("randomize", &AlignmentGenerator::GlobalShift::randomize)
        .def("__call__", &AlignmentGenerator::GlobalShift::operator());
  }

  {
    py::class_<AlignmentGenerator::GlobalRotation>(
        m, "AlignmentGeneratorGlobalRotation")
        .def(py::init<>())
        .def_readwrite("axis", &AlignmentGenerator::GlobalRotation::axis)
        .def_readwrite("angle", &AlignmentGenerator::GlobalRotation::angle)
        .def_readwrite("randomize",
                       &AlignmentGenerator::GlobalRotation::randomize)
        .def("__call__", &AlignmentGenerator::GlobalRotation::operator());
  }

  {
    py::class_<AlignmentGenerator::LocalRotation>(
        m, "AlignmentGeneratorLocalRotation")
        .def(py::init<>())
        .def_readwrite("axis", &AlignmentGenerator::LocalRotation::axis)
        .def_readwrite("angle", &AlignmentGenerator::LocalRotation::angle)
        .def_readwrite("randomize",
                       &AlignmentGenerator::LocalRotation::randomize)
        .def("__call__", &AlignmentGenerator::LocalRotation::operator());
  }

  {
    py::class_<AlignmentGenerator::LocalShift>(m,
                                               "AlignmentGeneratorLocalShift")
        .def(py::init<>())
        .def_readwrite("axisDirection",
                       &AlignmentGenerator::LocalShift::axisDirection)
        .def_readwrite("shift", &AlignmentGenerator::LocalShift::shift)
        .def_readwrite("randomize", &AlignmentGenerator::LocalShift::randomize)
        .def("__call__", &AlignmentGenerator::LocalShift::operator());
  }
}
