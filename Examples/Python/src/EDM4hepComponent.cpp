// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementInputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementOutputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMultiTrajectoryOutputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepParticleWriter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitWriter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimInputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepTrackInputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepTrackOutputConverter.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace Acts::Python;

PYBIND11_MODULE(ActsPythonBindingsEDM4hep, m) {
  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::EDM4hepSimInputConverter, m, "EDM4hepSimInputConverter",
      inputFrame, inputParticles, inputSimHits, outputParticlesGenerator,
      outputParticlesSimulation, outputSimHits, outputSimVertices,
      graphvizOutput, dd4hepDetector, trackingGeometry, sortSimHitsInTime);

  ACTS_PYTHON_DECLARE_WRITER(
      ActsExamples::EDM4hepSimHitWriter, m, "EDM4hepSimHitWriter", inputSimHits,
      inputParticles, outputPath, outputParticles, outputSimTrackerHits);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::EDM4hepMeasurementInputConverter,
                                m, "EDM4hepMeasurementInputConverter",
                                inputFrame, outputMeasurements,
                                outputMeasurementSimHitsMap, outputClusters);

  ACTS_PYTHON_DECLARE_ALGORITHM_CUSTOM(
      ActsExamples::EDM4hepMeasurementOutputConverter, m,
      "EDM4hepMeasurementOutputConverter", inputMeasurements, inputClusters,
      outputTrackerHitsPlane, outputTrackerHitsRaw) {
    alg.def_property_readonly(
        "collections",
        &ActsExamples::EDM4hepMeasurementOutputConverter::collections);
  };

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::EDM4hepParticleWriter, m,
                             "EDM4hepParticleWriter", inputParticles,
                             outputPath, outputParticles);

  ACTS_PYTHON_DECLARE_ALGORITHM_CUSTOM(
      ActsExamples::EDM4hepMultiTrajectoryOutputConverter, m,
      "EDM4hepMultiTrajectoryOutputConverter", inputTrajectories,
      inputMeasurementParticlesMap, outputTracks, Bz) {
    alg.def_property_readonly(
        "collections",
        &ActsExamples::EDM4hepMultiTrajectoryOutputConverter::collections);
  };

  ACTS_PYTHON_DECLARE_ALGORITHM_CUSTOM(
      ActsExamples::EDM4hepTrackOutputConverter, m,
      "EDM4hepTrackOutputConverter", inputTracks, outputTracks, Bz) {
    alg.def_property_readonly(
        "collections", &ActsExamples::EDM4hepTrackOutputConverter::collections);
  };

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::EDM4hepTrackInputConverter, m,
                                "EDM4hepTrackInputConverter", inputFrame,
                                inputTracks, outputTracks, Bz);
}
