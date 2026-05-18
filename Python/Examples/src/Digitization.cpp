// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/DigitizationConfigurator.hpp"
#include "ActsExamples/Digitization/DigitizationCoordinatesConverter.hpp"
#include "ActsExamples/Digitization/MuonSpacePointDigitizer.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <array>
#include <memory>
#include <tuple>
#include <utility>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace ActsPython {

void addDigitization(py::module& mex) {
  {
    auto [a, c] = declareAlgorithm<DigitizationAlgorithm, IAlgorithm>(
        mex, "DigitizationAlgorithm");

    ACTS_PYTHON_STRUCT(
        c, inputSimHits, outputMeasurements, outputClusters,
        outputMeasurementParticlesMap, outputMeasurementSimHitsMap,
        outputParticleMeasurementsMap, outputSimHitMeasurementsMap,
        surfaceByIdentifier, randomNumbers, doOutputCells, doClusterization,
        doMerge, mergeCommonCorner, minEnergyDeposit, digitizationConfigs,
        minMaxRetries);

    c.def_readonly("mergeNsigma", &DigitizationAlgorithm::Config::mergeNsigma);

    patchKwargsConstructor(c);

    auto cc = py::class_<DigiComponentsConfig>(mex, "DigiComponentsConfig")
                  .def(py::init<>());

    ACTS_PYTHON_STRUCT(cc, geometricDigiConfig, smearingDigiConfig);

    py::class_<DigiConfigContainer>(mex, "DigiConfigContainer")
        .def(py::init<std::vector<
                 std::pair<GeometryIdentifier, DigiComponentsConfig>>>());
  }

  {
    auto [alg, c] = declareAlgorithm<MuonSpacePointDigitizer, IAlgorithm>(
        mex, "MuonSpacePointDigitizer");

    ACTS_PYTHON_STRUCT(
        c, inputSimHits, inputParticles, outputSpacePoints, outputMeasurements,
        outputMeasurementParticlesMap, outputMeasurementSimHitsMap,
        outputParticleMeasurementsMap, outputSimHitMeasurementsMap,
        randomNumbers, trackingGeometry, digitizeTime, dumpVisualization,
        visualizationFunction, strawDeadTime, rpcDeadTime);

    ActsPython::patchKwargsConstructor(c);
  }

  {
    using DC = DigitizationConfigurator;
    auto dc = py::class_<DC>(mex, "DigitizationConfigurator").def(py::init<>());

    dc.def("__call__", &DC::operator());

    ACTS_PYTHON_STRUCT(dc, inputDigiComponents, compactify,
                       volumeLayerComponents, outputDigiComponents);
  }

  {
    py::class_<DigitizationCoordinatesConverter,
               std::shared_ptr<DigitizationCoordinatesConverter>>(
        mex, "DigitizationCoordinatesConverter")
        .def(py::init<DigitizationAlgorithm::Config&>(), py::arg("config"))
        .def_property_readonly("config",
                               &DigitizationCoordinatesConverter::config)
        .def("globalToLocal", &DigitizationCoordinatesConverter::globalToLocal)
        .def("localToGlobal", &DigitizationCoordinatesConverter::localToGlobal);
  }
}

}  // namespace ActsPython
