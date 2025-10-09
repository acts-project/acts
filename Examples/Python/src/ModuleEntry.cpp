// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/ActsVersion.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

#include <tuple>
#include <unordered_map>

#include <pybind11/detail/common.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pyerrors.h>

namespace py = pybind11;
using namespace ActsPython;

namespace ActsPython {
void addFramework(Context& ctx);

void addPropagation(Context& ctx);

void addAlignment(Context& ctx);

void addMaterialMapping(Context& ctx);
void addOutput(Context& ctx);
void addDetector(Context& ctx);
void addExampleAlgorithms(Context& ctx);
void addInput(Context& ctx);
void addGenerators(Context& ctx);
void addTruthTracking(Context& ctx);
void addTrackFitting(Context& ctx);
void addTrackFinding(Context& ctx);
void addTruthJet(Context& ctx);
void addVertexing(Context& ctx);
void addAmbiguityResolution(Context& ctx);
void addUtilities(Context& ctx);

void addRootInput(Context& ctx);
void addRootOutput(Context& ctx);

// Plugins
void addDigitization(Context& ctx);
void addPythia8(Context& ctx);
void addGeoModel(Context& ctx);
void addTGeo(Context& ctx);
void addJson(Context& ctx);
void addDetray(Context& ctx);
void addHepMC3(Context& ctx);
void addGnnTrackFinding(Context& ctx);
void addSvg(Context& ctx);
void addObj(Context& ctx);
void addOnnx(Context& ctx);
void addOnnxNeuralCalibrator(Context& ctx);
void addCovfie(Context& ctx);
void addTraccc(Context& ctx);
void addHashing(Context& ctx);

void addModuleEntry(Context& ctx) {
  addFramework(ctx);
  addOutput(ctx);

  addPropagation(ctx);
  addAlignment(ctx);
  addMaterialMapping(ctx);
  addDetector(ctx);
  addExampleAlgorithms(ctx);
  addInput(ctx);
  addGenerators(ctx);
  addTruthTracking(ctx);
  addTrackFitting(ctx);
  addTrackFinding(ctx);
  addTruthJet(ctx);
  addVertexing(ctx);
  addAmbiguityResolution(ctx);
  addUtilities(ctx);

  addDigitization(ctx);
  addPythia8(ctx);
  addJson(ctx);
  addGeoModel(ctx);
  addTGeo(ctx);
  addDetray(ctx);
  addHepMC3(ctx);
  addGnnTrackFinding(ctx);
  addObj(ctx);
  addSvg(ctx);
  addOnnx(ctx);
  addOnnxNeuralCalibrator(ctx);
  addCovfie(ctx);
  addTraccc(ctx);
  addHashing(ctx);

  addRootInput(ctx);
  addRootOutput(ctx);
}

}  // namespace ActsPython
