// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/DigitizationConfigurator.hpp"
#include "ActsExamples/Digitization/DigitizationCoordinatesConverter.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"

#include <array>
#include <memory>
#include <tuple>
#include <utility>

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
    using Config = ActsExamples::DigitizationAlgorithm::Config;

    auto a = py::class_<ActsExamples::DigitizationAlgorithm,
                        ActsExamples::IAlgorithm,
                        std::shared_ptr<ActsExamples::DigitizationAlgorithm>>(
                 mex, "DigitizationAlgorithm")
                 .def(py::init<Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def_property_readonly(
                     "config", &ActsExamples::DigitizationAlgorithm::config);

    auto c = py::class_<Config>(a, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(
        c, inputSimHits, outputMeasurements, outputClusters,
        outputMeasurementParticlesMap, outputMeasurementSimHitsMap,
        outputParticleMeasurementsMap, outputSimHitMeasurementsMap,
        surfaceByIdentifier, randomNumbers, doOutputCells, doClusterization,
        doMerge, minEnergyDeposit, digitizationConfigs);

    c.def_readonly("mergeNsigma", &Config::mergeNsigma);
    c.def_readonly("mergeCommonCorner", &Config::mergeCommonCorner);

    patchKwargsConstructor(c);

    auto cc = py::class_<DigiComponentsConfig>(mex, "DigiComponentsConfig")
                  .def(py::init<>());

    ACTS_PYTHON_STRUCT(cc, geometricDigiConfig, smearingDigiConfig);

    py::class_<DigiConfigContainer>(mex, "DigiConfigContainer")
        .def(py::init<std::vector<
                 std::pair<GeometryIdentifier, DigiComponentsConfig>>>());
  }

  {
    using DC = DigitizationConfigurator;
    auto dc = py::class_<DC>(mex, "DigitizationConfigurator").def(py::init<>());

    dc.def("__call__", &DC::operator());

    ACTS_PYTHON_STRUCT(dc, inputDigiComponents, compactify,
                       volumeLayerComponents, outputDigiComponents);
  }

  {
    py::class_<ActsExamples::DigitizationCoordinatesConverter,
               std::shared_ptr<ActsExamples::DigitizationCoordinatesConverter>>(
        mex, "DigitizationCoordinatesConverter")
        .def(py::init<ActsExamples::DigitizationAlgorithm::Config&>(),
             py::arg("config"))
        .def_property_readonly(
            "config", &ActsExamples::DigitizationCoordinatesConverter::config)
        .def("globalToLocal",
             &ActsExamples::DigitizationCoordinatesConverter::globalToLocal)
        .def("localToGlobal",
             &ActsExamples::DigitizationCoordinatesConverter::localToGlobal);
  }
}

}  // namespace Acts::Python
