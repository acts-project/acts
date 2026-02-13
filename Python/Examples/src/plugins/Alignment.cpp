// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"
#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"
#include "ActsExamples/DetectorCommons/AlignedTransformUpdater.hpp"
#include "ActsExamples/DetectorCommons/AlignmentDecorator.hpp"
#include "ActsExamples/DetectorCommons/AlignmentGenerator.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#ifdef ACTS_PLUGIN_DD4HEP
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#endif

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsAlignment, m) {
  // Bind AlignmentMask enum
  py::enum_<ActsAlignment::AlignmentMask>(m, "AlignmentMask",
                                          py::arithmetic())
      .value("None", ActsAlignment::AlignmentMask::None)
      .value("Center0", ActsAlignment::AlignmentMask::Center0)
      .value("Center1", ActsAlignment::AlignmentMask::Center1)
      .value("Center2", ActsAlignment::AlignmentMask::Center2)
      .value("Rotation0", ActsAlignment::AlignmentMask::Rotation0)
      .value("Rotation1", ActsAlignment::AlignmentMask::Rotation1)
      .value("Rotation2", ActsAlignment::AlignmentMask::Rotation2)
      .value("All", ActsAlignment::AlignmentMask::All)
      .export_values();

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
    py::class_<MutableGeoIdAlignmentStore, GeoIdAlignmentStore,
               std::shared_ptr<MutableGeoIdAlignmentStore>>(
        m, "MutableGeoIdAlignmentStore")
        .def(py::init<
             const std::unordered_map<GeometryIdentifier, Transform3>&>())
        .def("setTransform", &MutableGeoIdAlignmentStore::setTransform,
             py::arg("geoId"), py::arg("transform"),
             "Update or insert a transform for a given geometry ID")
        .def("getTransformMap",
             py::overload_cast<>(&MutableGeoIdAlignmentStore::getTransformMap),
             py::return_value_policy::reference_internal,
             "Get mutable access to the transform map");
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

  // Config
  py::class_<AA::Config>(m, "AlignmentAlgorithmConfig")
      .def(py::init<>())
      .def_readwrite("inputMeasurements", &AA::Config::inputMeasurements)
      .def_readwrite("inputProtoTracks", &AA::Config::inputProtoTracks)
      .def_readwrite("inputInitialTrackParameters",
                     &AA::Config::inputInitialTrackParameters)
      .def_readwrite("outputAlignmentParameters",
                     &AA::Config::outputAlignmentParameters)
      .def_readwrite("trackingGeometry", &AA::Config::trackingGeometry)
      .def_readwrite("alignedDetElements", &AA::Config::alignedDetElements)
      .def_readwrite("chi2ONdfCutOff", &AA::Config::chi2ONdfCutOff)
      .def_readwrite("align", &AA::Config::align)
      .def_readwrite("groups", &AA::Config::m_groups)
      .def_readwrite("alignedTransformUpdater",
                     &AA::Config::alignedTransformUpdater)
      .def_readwrite("deltaChi2ONdfCutOff", &AA::Config::deltaChi2ONdfCutOff)
      .def_readwrite("maxNumIterations", &AA::Config::maxNumIterations)
      .def_readwrite("maxNumTracks", &AA::Config::maxNumTracks)
      .def_property(
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
            cfg.iterationState.clear();
            for (auto item : d) {
              unsigned int iter = item.first.cast<unsigned int>();

              std::uint8_t mask_value;
              try {
                auto mask_enum =
                    item.second.cast<ActsAlignment::AlignmentMask>();
                mask_value = static_cast<std::uint8_t>(mask_enum);
              } catch (...) {
                mask_value = item.second.cast<std::uint8_t>();
              }

              cfg.iterationState[iter] = std::bitset<6>(mask_value);
            }
          });

  py::class_<AA::AlignmentFunction, std::shared_ptr<AA::AlignmentFunction>>(
      m, "AlignmentFunction");

  m.def("makeAlignmentFunction", &AA::makeAlignmentFunction,
        py::arg("trackingGeometry"), py::arg("magneticField"),
        py::arg("logLevel") = Acts::Logging::INFO);

  m.def("makeAlignedTransformUpdater",
        py::overload_cast<std::shared_ptr<MutableGeoIdAlignmentStore>>(
            &makeAlignedTransformUpdater),
        py::arg("store"),
        "Create an aligned transform updater with a mutable store");

  py::class_<AA, ActsExamples::IAlgorithm, std::shared_ptr<AA>>(
      m, "AlignmentAlgorithm")
      .def(py::init<AA::Config, Acts::Logging::Level>(), py::arg("config"),
           py::arg("level") = Acts::Logging::INFO);

  m.def(
      "surfacePlacement",
      [](const Acts::Surface& s) { return s.surfacePlacement(); },
      py::arg("surface"), py::return_value_policy::reference);
  
  // Deprecated: kept for backward compatibility
  m.def(
      "associatedDetectorElement",
      [](const Acts::Surface& s) {
        const auto* placement = s.surfacePlacement();
        return placement ? dynamic_cast<const Acts::DetectorElementBase*>(placement) : nullptr;
      },
      py::arg("surface"), py::return_value_policy::reference);
}
