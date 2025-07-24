// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/AmbiguityResolution/GreedyAmbiguityResolutionAlgorithm.hpp"
#include "ActsExamples/AmbiguityResolution/ScoreBasedAmbiguityResolutionAlgorithm.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;

namespace ActsPython {

/// This adds the ambiguity resolution algorithms to the examples module
/// @param mex the examples module
void addAmbiguityResolution(py::module_& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::GreedyAmbiguityResolutionAlgorithm, mex,
      "GreedyAmbiguityResolutionAlgorithm", inputTracks, outputTracks,
      maximumSharedHits, maximumIterations, nMeasurementsMin);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::ScoreBasedAmbiguityResolutionAlgorithm, mex,
      "ScoreBasedAmbiguityResolutionAlgorithm", inputTracks, configFile,
      outputTracks, minScore, minScoreSharedTracks, maxShared, minUnshared,
      maxSharedTracksPerMeasurement, useAmbiguityScoring);
}

}  // namespace ActsPython
