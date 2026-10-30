// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsAlignment/Kernel/Alignment.hpp"

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"
#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"
#include "ActsExamples/Alignment/ExampleAlignedTransformUpdater.hpp"
#include "ActsExamples/DetectorCommons/AlignmentDecorator.hpp"
#include "ActsExamples/DetectorCommons/AlignmentGenerator.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsAlignment, m) {
  py::enum_<ActsAlignment::AlignmentMask>(m, "AlignmentMask", py::arithmetic())
      .value("eNone", ActsAlignment::AlignmentMask::None)
      .value("Center0", ActsAlignment::AlignmentMask::Center0)
      .value("Center1", ActsAlignment::AlignmentMask::Center1)
      .value("Center2", ActsAlignment::AlignmentMask::Center2)
      .value("Rotation0", ActsAlignment::AlignmentMask::Rotation0)
      .value("Rotation1", ActsAlignment::AlignmentMask::Rotation1)
      .value("Rotation2", ActsAlignment::AlignmentMask::Rotation2)
      .value("All", ActsAlignment::AlignmentMask::All);

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
        cfg.iterationState.clear();
        for (auto item : d) {
          unsigned int iter = item.first.cast<unsigned int>();

          std::uint8_t mask_value;
          if (py::isinstance<ActsAlignment::AlignmentMask>(item.second)) {
            mask_value = static_cast<std::uint8_t>(
                item.second.cast<ActsAlignment::AlignmentMask>());
          } else if (py::isinstance<py::int_>(item.second)) {
            mask_value = item.second.cast<std::uint8_t>();
          } else {
            throw py::type_error(
                "iterationState value must be AlignmentMask or int");
          }

          cfg.iterationState[iter] = std::bitset<6>(mask_value);
        }
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
