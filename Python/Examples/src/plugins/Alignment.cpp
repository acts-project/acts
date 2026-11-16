// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Alignment.hpp"

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Geometry/AlignableStructure.hpp"
#include "ActsAlignment/Geometry/AlignmentHierarchy.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"
#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"
#include "ActsExamples/Alignment/ExampleAlignedTransformUpdater.hpp"
#include "ActsExamples/DetectorCommons/AlignmentDecorator.hpp"
#include "ActsExamples/DetectorCommons/AlignmentGenerator.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <bitset>
#include <map>

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

    py::enum_<AlignmentDecorator::Target>(ad, "Target")
        .value("eSim", AlignmentDecorator::Target::eSim)
        .value("eReco", AlignmentDecorator::Target::eReco)
        .value("eBoth", AlignmentDecorator::Target::eBoth);

    auto c =
        py::class_<AlignmentDecorator::Config>(ad, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, target, iovStores, nominalStore, garbageCollection,
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
             const std::unordered_map<GeometryIdentifier, Transform3>&>())
        .def("setTransform", &GeoIdAlignmentStore::setTransform,
             py::arg("geoId"), py::arg("transform"),
             "Update or insert a transform for a given geometry ID")
        .def("getTransformMap", &GeoIdAlignmentStore::getTransformMap,
             py::return_value_policy::reference_internal,
             "Return a new dict whose Transform3 values are live views of the "
             "stored transforms");
  }

  {
    py::class_<AlignableStructure, std::shared_ptr<AlignableStructure>>(
        m, "AlignableStructure")
        .def(py::init<Acts::GeometryIdentifier>(), py::arg("id"))
        .def(py::init<Acts::GeometryIdentifier, const Acts::Transform3&>(),
             py::arg("id"), py::arg("transform"))
        .def(py::init<Acts::GeometryIdentifier,
                      AlignableStructure::TransformProvider>(),
             py::arg("id"), py::arg("transformProvider"))
        .def("setTransform",
             py::overload_cast<const Acts::Transform3&>(
                 &AlignableStructure::setTransform),
             py::arg("transform"))
        .def("setTransform",
             py::overload_cast<AlignableStructure::TransformProvider>(
                 &AlignableStructure::setTransform),
             py::arg("transformProvider"))
        .def_property_readonly("hasTransform",
                               &AlignableStructure::hasTransform)
        .def("localToGlobalTransform",
             &AlignableStructure::localToGlobalTransform, py::arg("gctx"))
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
            [](const AlignableStructure& s) {
              const auto children = s.children();
              return std::vector<std::shared_ptr<AlignableStructure>>(
                  children.begin(), children.end());
            })
        .def_property(
            "alignmentMask",
            [](const AlignableStructure& s) { return s.alignmentMask(); },
            [](AlignableStructure& s, AlignmentMask mask) {
              s.alignmentMask() = mask;
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
                 static_cast<const AlignableStructure* (AlignmentHierarchy::*)(
                     const Acts::SurfacePlacementBase&) const>(
                     &AlignmentHierarchy::structureFor),
                 py::arg("detElement"),
                 py::return_value_policy::reference_internal)
            .def("structureFor",
                 static_cast<const AlignableStructure* (AlignmentHierarchy::*)(
                     const Acts::Surface&) const>(
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

  using AA = ActsExamples::AlignmentAlgorithm;

  auto c =
      py::class_<AA::Config>(m, "AlignmentAlgorithmConfig").def(py::init<>());

  ACTS_PYTHON_STRUCT(c, inputMeasurements, inputProtoTracks,
                     inputInitialTrackParameters, outputAlignmentParameters,
                     alignedDetElements, chi2ONdfCutOff, align,
                     alignedTransformUpdater, deltaChi2ONdfCutOff,
                     maxNumIterations, maxNumTracks);

  c.def_property(
      "iterationState",
      [](const AA::Config& cfg) {
        py::dict result;
        for (const auto& [iter, mask] : cfg.iterationState) {
          result[py::cast(iter)] =
              py::cast(static_cast<std::uint8_t>(mask.to_ulong()));
        }
        return result;
      },
      [](AA::Config& cfg, py::dict d) {
        // Convert first. A bad key or value throws before the stored map is
        // replaced, so a rejected assignment keeps the previous configuration.
        std::map<unsigned int, std::bitset<6>> parsed;
        for (auto item : d) {
          unsigned int iter = 0;
          try {
            iter = item.first.cast<unsigned int>();
          } catch (const py::cast_error&) {
            throw py::type_error("iterationState key must be an unsigned int");
          }

          std::uint8_t mask_value = 0;
          if (py::isinstance<ActsAlignment::AlignmentMask>(item.second)) {
            mask_value = static_cast<std::uint8_t>(
                item.second.cast<ActsAlignment::AlignmentMask>());
          } else if (py::isinstance<py::int_>(item.second)) {
            mask_value = item.second.cast<std::uint8_t>();
          } else {
            throw py::type_error(
                "iterationState value must be AlignmentMask or int");
          }

          parsed.emplace(iter, std::bitset<6>(mask_value));
        }
        cfg.iterationState = std::move(parsed);
      });

  py::class_<AA::AlignmentFunction, std::shared_ptr<AA::AlignmentFunction>>(
      m, "AlignmentFunction");

  m.def("makeAlignmentFunction", &AA::makeAlignmentFunction,
        py::arg("trackingGeometry"), py::arg("magneticField"));

  m.def("makeAlignedTransformUpdater", &makeAlignedTransformUpdater,
        py::arg("store"),
        "Create an aligned transform updater with an alignment store");

  py::class_<AA, ActsExamples::IAlgorithm, std::shared_ptr<AA>>(
      m, "AlignmentAlgorithm")
      .def(py::init([](AA::Config cfg, Acts::Logging::Level level) {
             return std::make_shared<AA>(
                 std::move(cfg), getDefaultLogger("AlignmentAlgorithm", level));
           }),
           py::arg("config"), py::arg("level") = Acts::Logging::INFO);
}
