// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/ActsVersion.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Plugins/FpeMonitoring/FpeMonitor.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <pybind11/detail/common.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pyerrors.h>

namespace py = pybind11;
using namespace Acts::Python;

namespace Acts::Python {
void addUnits(Context& ctx);
void addFramework(Context& ctx);
void addLogging(Context& ctx);
void addPdgParticle(Context& ctx);
void addAlgebra(Context& ctx);
void addBinning(Context& ctx);
void addEventData(Context& ctx);

void addPropagation(Context& ctx);

void addGeometry(Context& ctx);
void addExperimentalGeometry(Context& ctx);

void addMagneticField(Context& ctx);

void addMaterial(Context& ctx);
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

// Plugins
void addDigitization(Context& ctx);
void addPythia8(Context& ctx);
void addJson(Context& ctx);
void addHepMC3(Context& ctx);
void addExaTrkXTrackFinding(Context& ctx);
void addEDM4hep(Context& ctx);
void addSvg(Context& ctx);
void addObj(Context& ctx);
void addOnnx(Context& ctx);
void addOnnxMlpack(Context& ctx);
void addOnnxNeuralCalibrator(Context& ctx);

}  // namespace Acts::Python

PYBIND11_MODULE(ActsPythonBindings, m) {
  Acts::Python::Context ctx;
  ctx.modules["main"] = m;
  auto mex = m.def_submodule("_examples");
  ctx.modules["examples"] = mex;
  auto prop = m.def_submodule("_propagator");
  ctx.modules["propagation"] = prop;
  m.doc() = "Acts";

  m.attr("__version__") =
      std::tuple{Acts::VersionMajor, Acts::VersionMinor, Acts::VersionPatch};

  addUnits(ctx);
  addFramework(ctx);
  addLogging(ctx);
  addPdgParticle(ctx);
  addAlgebra(ctx);
  addBinning(ctx);
  addEventData(ctx);

  addPropagation(ctx);
  addGeometry(ctx);
  addExperimentalGeometry(ctx);

  addMagneticField(ctx);
  addMaterial(ctx);
  addOutput(ctx);
  addDetector(ctx);
  addExampleAlgorithms(ctx);
  addInput(ctx);
  addGenerators(ctx);
  addTruthTracking(ctx);
  addTrackFitting(ctx);
  addTrackFinding(ctx);
  addVertexing(ctx);
  addAmbiguityResolution(ctx);

  addDigitization(ctx);
  addPythia8(ctx);
  addJson(ctx);
  addHepMC3(ctx);
  addExaTrkXTrackFinding(ctx);
  addEDM4hep(ctx);
  addObj(ctx);
  addSvg(ctx);
  addOnnx(ctx);
  addOnnxMlpack(ctx);
  addOnnxNeuralCalibrator(ctx);
}
