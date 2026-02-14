// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementInputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementOutputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMultiTrajectoryOutputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepParticleOutputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitOutputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimInputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepTrackInputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepTrackOutputConverter.hpp"
#include "ActsExamples/Io/Podio/PodioMeasurementInputConverter.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"
#include "ActsExamples/Io/Podio/PodioReader.hpp"
#include "ActsExamples/Io/Podio/PodioWriter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <podio/CollectionBase.h>
#include <podio/Frame.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsPython;
using namespace ActsExamples;

PYBIND11_MODULE(ActsExamplesPythonBindingsEDM4hep, m) {
  ACTS_PYTHON_DECLARE_READER(PodioReader, m, "PodioReader", inputPath,
                             outputFrame, category);

  ACTS_PYTHON_DECLARE_WRITER(PodioWriter, m, "PodioWriter", inputFrame,
                             outputPath, category, collections,
                             separateFilesPerThread);

  py::class_<PodioOutputConverter, IAlgorithm,
             std::shared_ptr<PodioOutputConverter>>(m, "PodioOutputConverter")
      .def_property_readonly("collections", &PodioOutputConverter::collections);

  py::class_<PodioInputConverter, IAlgorithm,
             std::shared_ptr<PodioInputConverter>>(m, "PodioInputConverter");

  {
    auto [alg, config] =
        declareAlgorithm<PodioMeasurementInputConverter, PodioInputConverter>(
            m, "PodioMeasurementInputConverter");
    ACTS_PYTHON_STRUCT(
        config, inputMeasurements, inputFrame, outputMeasurements,
        outputMeasurementParticlesMap, outputMeasurementSimHitsMap,
        outputParticleMeasurementsMap, outputSimHitMeasurementsMap,
        inputSimHits, inputSimHitAssociation);
  }

  {
    auto [alg, config] = declareAlgorithm<EDM4hepSimInputConverter, IAlgorithm>(
        m, "EDM4hepSimInputConverter");
    ACTS_PYTHON_STRUCT(
        config, inputFrame, inputParticles, inputSimHits,
        outputParticlesGenerator, outputParticlesSimulation, outputSimHits,
        outputSimHitAssociation, outputSimVertices, dd4hepDetector,
        trackingGeometry, sortSimHitsInTime, particleRMin, particleRMax,
        particleZMin, particleZMax, particlePtMin, particlePtMax);

    using Config = EDM4hepSimInputConverter::Config;
    pythonRangeProperty(config, "particleR", &Config::particleRMin,
                        &Config::particleRMax);

    pythonRangeProperty(config, "particleZ", &Config::particleZMin,
                        &Config::particleZMax);

    pythonRangeProperty(config, "pt", &Config::particlePtMin,
                        &Config::particlePtMax);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(EDM4hepTrackInputConverter, m,
                                "EDM4hepTrackInputConverter", inputFrame,
                                inputTracks, outputTracks, Bz);

  {
    auto [alg, config] =
        declareAlgorithm<EDM4hepSimHitOutputConverter, PodioOutputConverter>(
            m, "EDM4hepSimHitOutputConverter");
    ACTS_PYTHON_STRUCT(config, inputSimHits, inputParticles, outputParticles,
                       outputSimTrackerHits);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(EDM4hepMeasurementInputConverter, m,
                                "EDM4hepMeasurementInputConverter", inputFrame,
                                outputMeasurements, outputMeasurementSimHitsMap,
                                outputClusters);

  {
    auto [alg, config] = declareAlgorithm<EDM4hepMeasurementOutputConverter,
                                          PodioOutputConverter>(
        m, "EDM4hepMeasurementOutputConverter");
    ACTS_PYTHON_STRUCT(config, inputMeasurements, inputClusters,
                       outputTrackerHitsPlane, outputTrackerHitsRaw);
  }

  {
    auto [alg, config] =
        declareAlgorithm<EDM4hepParticleOutputConverter, PodioOutputConverter>(
            m, "EDM4hepParticleOutputConverter");
    ACTS_PYTHON_STRUCT(config, inputParticles, outputParticles);
  }

  {
    auto [alg, config] = declareAlgorithm<EDM4hepMultiTrajectoryOutputConverter,
                                          PodioOutputConverter>(
        m, "EDM4hepMultiTrajectoryOutputConverter");
    ACTS_PYTHON_STRUCT(config, inputTrajectories, inputMeasurementParticlesMap,
                       outputTracks, Bz);
  }

  {
    auto [alg, config] =
        declareAlgorithm<EDM4hepTrackOutputConverter, PodioOutputConverter>(
            m, "EDM4hepTrackOutputConverter");
    ACTS_PYTHON_STRUCT(config, inputTracks, outputTracks, Bz);
  }
}
