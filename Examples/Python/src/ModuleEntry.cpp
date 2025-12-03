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
void addVertexing(Context& ctx);
void addAmbiguityResolution(Context& ctx);
void addUtilities(Context& ctx);

// Plugins
void addDigitization(Context& ctx);
void addObj(Context& ctx);

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
  addVertexing(ctx);
  addAmbiguityResolution(ctx);
  addUtilities(ctx);

  addDigitization(ctx);
  addObj(ctx);
}

}  // namespace ActsPython
