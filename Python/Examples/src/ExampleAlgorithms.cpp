// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/Fatras/FatrasSimulation.hpp"
#include "ActsExamples/Printers/ParticlesPrinter.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Utilities/TrackSelectorAlgorithm.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace ActsPython {

void addExampleAlgorithms(py::module& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::FatrasSimulation, mex, "FatrasSimulation", inputParticles,
      outputParticles, outputSimHits, randomNumbers, trackingGeometry,
      magneticField, pMin, emScattering, emEnergyLossIonisation,
      emEnergyLossRadiation, emPhotonConversion, generateHitsOnSensitive,
      generateHitsOnMaterial, generateHitsOnPassive, averageHitsPerParticle);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::ParticlesPrinter, mex,
                                "ParticlesPrinter", inputParticles);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackParametersPrinter, mex,
                                "TrackParametersPrinter", inputTrackParameters);

  {
    using Alg = ActsExamples::TrackSelectorAlgorithm;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, IAlgorithm, std::shared_ptr<Alg>>(
                   mex, "TrackSelectorAlgorithm")
                   .def(py::init<const Alg::Config&, Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, inputTracks, outputTracks, selectorConfig);
  }
}
}  // namespace ActsPython
