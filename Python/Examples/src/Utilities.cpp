// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/MeasurementMapSelector.hpp"
#include "ActsExamples/Utilities/ProtoTracksToSeeds.hpp"
#include "ActsExamples/Utilities/ProtoTracksToTracks.hpp"
#include "ActsExamples/Utilities/SeedsToProtoTracks.hpp"
#include "ActsExamples/Utilities/TracksToParameters.hpp"
#include "ActsExamples/Utilities/TracksToTrajectories.hpp"
#include "ActsExamples/Utilities/TrajectoriesToProtoTracks.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addUtilities(py::module& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(TrajectoriesToProtoTracks, mex,
                                "TrajectoriesToProtoTracks", inputTrajectories,
                                outputProtoTracks);

  ACTS_PYTHON_DECLARE_ALGORITHM(TracksToTrajectories, mex,
                                "TracksToTrajectories", inputTracks,
                                outputTrajectories);

  ACTS_PYTHON_DECLARE_ALGORITHM(TracksToParameters, mex, "TracksToParameters",
                                inputTracks, outputTrackParameters);

  ACTS_PYTHON_DECLARE_ALGORITHM(SeedsToProtoTracks, mex, "SeedsToProtoTracks",
                                inputSeeds, outputProtoTracks);

  ACTS_PYTHON_DECLARE_ALGORITHM(ProtoTracksToSeeds, mex, "ProtoTracksToSeeds",
                                inputProtoTracks, inputSpacePoints, outputSeeds,
                                outputProtoTracks);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      MeasurementMapSelector, mex, "MeasurementMapSelector", inputMeasurements,
      inputMeasurementParticleMap, outputMeasurementParticleMap,
      geometrySelection);

  ACTS_PYTHON_DECLARE_ALGORITHM(ProtoTracksToTracks, mex, "ProtoTracksToTracks",
                                inputMeasurements, inputProtoTracks,
                                inputTrackParameters, outputTracks);
}

}  // namespace ActsPython
