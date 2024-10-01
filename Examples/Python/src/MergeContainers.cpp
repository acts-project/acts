// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/MergeContainers/MergeContainersAlgorithm.hpp"
#include "ActsExamples/Matching/MatchingAlgorithm.hpp"
#include "ActsExamples/TransportParticles/TransportParticlesAlgorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;

namespace Acts::Python {

void addMergeContainers(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::MergeContainersAlgorithm, mex,
      "MergeContainersAlgorithm", inputTrackParameters, outputTrackParameters, inputTracks, outputTracks);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::MatchingAlgorithm, mex, "MatchingAlgorithm",
      inputTrackParametersMS,
      inputTrackParametersVT,
      inputTrackContainerMS,
      inputTrackContainerVT,
      outputTrackParameters,
      outputTracks,
      outputMatchedTracks,
      inputParticles,
      inputMeasurementParticlesMapVT,
      inputMeasurementParticlesMapMS,
      useRecVtx, px, py, pz, chi2max, fit,
      trackingGeometry, magneticField);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TransportParticles, mex, "TransportParticlesAlgorithm",
      inputTrackParameters,
      inputTrackContainer,
      outputTracks,
      inputVertices,
      useRecVtx, px, py, pz,
      trackingGeometry, magneticField);

}

}  // namespace Acts::Python
