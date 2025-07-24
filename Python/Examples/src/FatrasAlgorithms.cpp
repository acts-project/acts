// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Fatras/FatrasSimulation.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include "ActsPython/Utilities/Patchers.hpp"

#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace ActsPython {

/// This adds the Fatras algorithms to the examples module
/// @param mex the examples module
void addFatrasAlgorithms(py::module_& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::FatrasSimulation, mex, "FatrasSimulation", inputParticles,
      outputParticles, outputSimHits, randomNumbers, trackingGeometry,
      magneticField, pMin, emScattering, emEnergyLossIonisation,
      emEnergyLossRadiation, emPhotonConversion, generateHitsOnSensitive,
      generateHitsOnMaterial, generateHitsOnPassive, averageHitsPerParticle);
}
}  // namespace ActsPython
