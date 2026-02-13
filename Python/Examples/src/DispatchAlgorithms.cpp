// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DispatchAlgorithms/DispatchEdm.hpp"
#include "ActsExamples/DispatchAlgorithms/PatternDispatchAlgorithm.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addDispatchAlgorithms(py::module& mex) {
  // Dispatch event data model as structurs

  auto dmeas =
      py::class_<DispatchMeasurements, std::shared_ptr<DispatchMeasurements>>(
          mex, "DispatchMeasurements")
          .def(py::init<>());
  ACTS_PYTHON_STRUCT(dmeas, clusterIndices, clusterGeoIds, x, y, z, lx, ly,
                     covLxx, covLyy);

  auto dpart =
      py::class_<DispatchParticles, std::shared_ptr<DispatchParticles>>(
          mex, "DispatchParticles")
          .def(py::init<>());
  ACTS_PYTHON_STRUCT(dpart, particlePdgs, px, py, pz, vx, vy, vz);

  auto dtrack =
      py::class_<DispatchTrack, std::shared_ptr<DispatchTrack>>(
          mex, "DispatchTrack")
          .def(py::init<>());
  ACTS_PYTHON_STRUCT(dtrack, parameters, covariances, measurementIndices);

  auto dmap = py::class_<DispatchParticleMeasurementsMap,
                         std::shared_ptr<DispatchParticleMeasurementsMap>>(
                  mex, "DispatchParticleMeasurementsMap")
                  .def(py::init<>());

  // Pattern recognition algorithm
  {
    ACTS_PYTHON_DECLARE_ALGORITHM(
        ActsExamples::PatternDispatchAlgorithm, mex, "PatternDispatchAlgorithm",
        patternFunction, trackingGeometry, inputParticles,
        inputParticleMeasurementsMap, inputMeasurements, outputProtoTracks);
  }
}
}  // namespace ActsPython
