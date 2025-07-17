// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPython/Module/Entries.hpp"
#include "ActsPython/Utilities/Context.hpp"

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

namespace ActsPython {
void addFramework(Context& ctx);

void addPropagation(Context& ctx);
void addMagneticFieldLegacy(Context& ctx);

void addMaterialLegacy(Context& ctx);
void addOutput(Context& ctx);
void addDetector(Context& ctx);
void addExampleAlgorithms(Context& ctx);
void addInput(Context& ctx);
void addGenerators(Context& ctx);
void addTruthTracking(Context& ctx);
void addTrackFitting(Context& ctx);
void addTrackFindingLegacy(Context& ctx);
void addTruthJet(Context& ctx);
void addVertexingLegacy(Context& ctx);
void addAmbiguityResolution(Context& ctx);
void addUtilitiesLegacy(Context& ctx);

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
void addExaTrkXTrackFinding(Context& ctx);
void addSvg(Context& ctx);
void addObj(Context& ctx);
void addOnnx(Context& ctx);
void addOnnxNeuralCalibrator(Context& ctx);
void addCovfie(Context& ctx);
void addTraccc(Context& ctx);
void addHashing(Context& ctx);

}  // namespace ActsPython

void ActsPython::addLegacyExamplesModule(Context& ctx) {
  auto& m = ctx.get("main");

  auto mex = m.def_submodule("_examples");
  ctx.modules["examples"] = mex;
  auto prop = m.def_submodule("_propagator");
  ctx.modules["propagation"] = prop;

  addFramework(ctx);
  addOutput(ctx);

  addPropagation(ctx);

  addMagneticFieldLegacy(ctx);
  addMaterialLegacy(ctx);
  addDetector(ctx);
  addExampleAlgorithms(ctx);
  addInput(ctx);
  addGenerators(ctx);
  addTruthTracking(ctx);
  addTrackFitting(ctx);
  addTrackFindingLegacy(ctx);
  addTruthJet(ctx);
  addVertexingLegacy(ctx);
  addAmbiguityResolution(ctx);
  addUtilitiesLegacy(ctx);

  addDigitization(ctx);
  addPythia8(ctx);
  addJson(ctx);
  addGeoModel(ctx);
  addTGeo(ctx);
  addDetray(ctx);
  addHepMC3(ctx);
  addExaTrkXTrackFinding(ctx);
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
