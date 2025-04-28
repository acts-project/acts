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
#include "ActsExamples/Io/EDM4hep/EDM4hepOutputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepParticleOutputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitOutputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimInputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepTrackInputConverter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepTrackOutputConverter.hpp"

#include <podio/CollectionBase.h>
#include <podio/Frame.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace Acts::Python;
using namespace ActsExamples;

template <typename A, typename B = IAlgorithm>
auto declareAlgorithm(py::module_& m, const char* name) {
  using Config = typename A::Config;
  auto alg = py::class_<A, B, std::shared_ptr<A>>(m, name)
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def_property_readonly("config", &A::config);
  auto c = py::class_<Config>(alg, "Config").def(py::init<>());
  return std::tuple{alg, c};
}

PYBIND11_MODULE(ActsPythonBindingsEDM4hep, m) {
  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::EDM4hepSimInputConverter, m, "EDM4hepSimInputConverter",
      inputFrame, inputParticles, inputSimHits, outputParticlesGenerator,
      outputParticlesSimulation, outputSimHits, outputSimVertices,
      graphvizOutput, dd4hepDetector, trackingGeometry, sortSimHitsInTime);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::EDM4hepTrackInputConverter, m,
                                "EDM4hepTrackInputConverter", inputFrame,
                                inputTracks, outputTracks, Bz);

  py::class_<EDM4hepOutputConverter, IAlgorithm,
             std::shared_ptr<EDM4hepOutputConverter>>(m,
                                                      "EDM4hepOutputConverter")
      .def_property_readonly("collections",
                             &EDM4hepOutputConverter::collections);

  {
    auto [alg, config] =
        declareAlgorithm<EDM4hepSimHitOutputConverter, EDM4hepOutputConverter>(
            m, "EDM4hepSimHitOutputConverter");
    ACTS_PYTHON_STRUCT(config, inputSimHits, inputParticles, outputParticles,
                       outputSimTrackerHits);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::EDM4hepMeasurementInputConverter,
                                m, "EDM4hepMeasurementInputConverter",
                                inputFrame, outputMeasurements,
                                outputMeasurementSimHitsMap, outputClusters);

  {
    auto [alg, config] = declareAlgorithm<EDM4hepMeasurementOutputConverter,
                                          EDM4hepOutputConverter>(
        m, "EDM4hepMeasurementOutputConverter");
    ACTS_PYTHON_STRUCT(config, inputMeasurements, inputClusters,
                       outputTrackerHitsPlane, outputTrackerHitsRaw);
  }

  {
    auto [alg, config] = declareAlgorithm<EDM4hepParticleOutputConverter,
                                          EDM4hepOutputConverter>(
        m, "EDM4hepParticleOutputConverter");
    ACTS_PYTHON_STRUCT(config, inputParticles, outputParticles);
  }

  {
    auto [alg, config] = declareAlgorithm<EDM4hepMultiTrajectoryOutputConverter,
                                          EDM4hepOutputConverter>(
        m, "EDM4hepMultiTrajectoryOutputConverter");
    ACTS_PYTHON_STRUCT(config, inputTrajectories, inputMeasurementParticlesMap,
                       outputTracks, Bz);
  }

  {
    auto [alg, config] =
        declareAlgorithm<EDM4hepTrackOutputConverter, EDM4hepOutputConverter>(
            m, "EDM4hepTrackOutputConverter");
    ACTS_PYTHON_STRUCT(config, inputTracks, outputTracks, Bz);
  }
}
