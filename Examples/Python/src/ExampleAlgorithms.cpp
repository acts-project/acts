// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Fatras/FatrasSimulation.hpp"
#include "ActsExamples/Io/Json/JsonGeometryList.hpp"
#include "ActsExamples/Printers/HitsPrinter.hpp"
#include "ActsExamples/Printers/ParticlesPrinter.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addExampleAlgorithms(Context& ctx) {
  auto mex = ctx.get("examples");

  mex.def("readJsonGeometryList", ActsExamples::readJsonGeometryList);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::FatrasSimulation, mex, "FatrasSimulation", inputParticles,
      outputParticlesInitial, outputParticlesFinal, outputSimHits,
      imputParametrisationNuclearInteraction, randomNumbers, trackingGeometry,
      magneticField, pMin, emScattering, emEnergyLossIonisation,
      emEnergyLossRadiation, emPhotonConversion, generateHitsOnSensitive,
      generateHitsOnMaterial, generateHitsOnPassive, averageHitsPerParticle);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::ParticlesPrinter, mex,
                                "ParticlesPrinter", inputParticles);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::HitsPrinter, mex, "HitsPrinter", inputClusters,
      inputMeasurementParticlesMap, inputHitIds, selectIndexStart,
      selectIndexLength, selectVolume, selectLayer, selectModule);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackParametersPrinter, mex,
                                "TrackParametersPrinter", inputTrackParameters);
}

}  // namespace Acts::Python
