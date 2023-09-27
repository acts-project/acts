// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Digitization/PlanarModuleStepper.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/DigitizationConfigurator.hpp"
#include "ActsExamples/Digitization/PlanarSteppingAlgorithm.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addDigitization(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  mex.def("readDigiConfigFromJson", ActsExamples::readDigiConfigFromJson);
  mex.def("writeDigiConfigToJson", ActsExamples::writeDigiConfigToJson);

  {
    using Config = ActsExamples::DigitizationConfig;

    py::class_<ActsExamples::DigitizationAlgorithm, ActsExamples::IAlgorithm,
               std::shared_ptr<ActsExamples::DigitizationAlgorithm>>(
        mex, "DigitizationAlgorithm")
        .def(py::init<Config&, Acts::Logging::Level>(), py::arg("config"),
             py::arg("level"))
        .def_property_readonly("config",
                               &ActsExamples::DigitizationAlgorithm::config);

    auto c = py::class_<Config>(mex, "DigitizationConfig")
                 .def(py::init<Acts::GeometryHierarchyMap<
                          ActsExamples::DigiComponentsConfig>>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(outputSourceLinks);
    ACTS_PYTHON_MEMBER(outputMeasurements);
    ACTS_PYTHON_MEMBER(outputClusters);
    ACTS_PYTHON_MEMBER(outputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(outputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(randomNumbers);
    ACTS_PYTHON_MEMBER(doMerge);
    ACTS_PYTHON_MEMBER(digitizationConfigs);
    ACTS_PYTHON_STRUCT_END();

    c.def_readonly("mergeNsigma", &Config::mergeNsigma);
    c.def_readonly("mergeCommonCorner", &Config::mergeCommonCorner);

    patchKwargsConstructor(c);

    py::class_<DigiComponentsConfig>(mex, "DigiComponentsConfig");

    py::class_<Acts::GeometryHierarchyMap<ActsExamples::DigiComponentsConfig>>(
        mex, "GeometryHierarchyMap_DigiComponentsConfig")
        .def(py::init<std::vector<
                 std::pair<GeometryIdentifier, DigiComponentsConfig>>>());
  }

  {
    using Alg = ActsExamples::PlanarSteppingAlgorithm;

    auto alg = py::class_<Alg, ActsExamples::IAlgorithm, std::shared_ptr<Alg>>(
                   mex, "PlanarSteppingAlgorithm")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Alg::Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Alg::Config);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(outputClusters);
    ACTS_PYTHON_MEMBER(outputSourceLinks);
    ACTS_PYTHON_MEMBER(outputDigiSourceLinks);
    ACTS_PYTHON_MEMBER(outputMeasurements);
    ACTS_PYTHON_MEMBER(outputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(outputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(planarModuleStepper);
    ACTS_PYTHON_MEMBER(randomNumbers);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    py::class_<PlanarModuleStepper, std::shared_ptr<PlanarModuleStepper>>(
        m, "PlanarModuleStepper")
        .def(py::init<>())
        .def(py::init([](Logging::Level level) {
          return std::make_shared<PlanarModuleStepper>(
              getDefaultLogger("PlanarModuleStepper", level));
        }));
  }

  {
    using DC = DigitizationConfigurator;
    auto dc = py::class_<DC>(mex, "DigitizationConfigurator").def(py::init<>());

    dc.def("__call__", &DC::operator());

    ACTS_PYTHON_STRUCT_BEGIN(dc, DC);
    ACTS_PYTHON_MEMBER(inputDigiComponents);
    ACTS_PYTHON_MEMBER(compactify);
    ACTS_PYTHON_MEMBER(volumeLayerComponents);
    ACTS_PYTHON_MEMBER(outputDigiComponents);
    ACTS_PYTHON_STRUCT_END();
  }
}

}  // namespace Acts::Python
