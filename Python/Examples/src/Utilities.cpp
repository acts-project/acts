// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/MeasurementMapSelector.hpp"
#include "ActsExamples/Utilities/PrototracksToParameters.hpp"
#include "ActsExamples/Utilities/PrototracksToSeeds.hpp"
#include "ActsExamples/Utilities/PrototracksToTracks.hpp"
#include "ActsExamples/Utilities/SeedsToPrototracks.hpp"
#include "ActsExamples/Utilities/TracksToParameters.hpp"
#include "ActsExamples/Utilities/TracksToTrajectories.hpp"
#include "ActsExamples/Utilities/TrajectoriesToPrototracks.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addUtilities(py::module& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(TrajectoriesToPrototracks, mex,
                                "TrajectoriesToPrototracks", inputTrajectories,
                                outputProtoTracks);

  ACTS_PYTHON_DECLARE_ALGORITHM(TracksToTrajectories, mex,
                                "TracksToTrajectories", inputTracks,
                                outputTrajectories);

  ACTS_PYTHON_DECLARE_ALGORITHM(TracksToParameters, mex, "TracksToParameters",
                                inputTracks, outputTrackParameters);

  ACTS_PYTHON_DECLARE_ALGORITHM(SeedsToPrototracks, mex, "SeedsToPrototracks",
                                inputSeeds, outputProtoTracks);

  ACTS_PYTHON_DECLARE_ALGORITHM(PrototracksToSeeds, mex, "PrototracksToSeeds",
                                inputProtoTracks, inputSpacePoints, outputSeeds,
                                outputProtoTracks);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      MeasurementMapSelector, mex, "MeasurementMapSelector", inputMeasurements,
      inputMeasurementParticleMap, outputMeasurementParticleMap,
      geometrySelection);

  ACTS_PYTHON_DECLARE_ALGORITHM(PrototracksToTracks, mex, "PrototracksToTracks",
                                inputMeasurements, inputProtoTracks,
                                inputTrackParameters, outputTracks);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      PrototracksToParameters, mex, "PrototracksToParameters", inputProtoTracks,
      inputSpacePoints, outputSeeds, outputParameters, outputProtoTracks,
      geometry, magneticField, buildTightSeeds);
}

}  // namespace ActsPython
