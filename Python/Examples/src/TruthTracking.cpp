// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/TruthTracking/HitSelector.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/TruthTracking/ParticleTrackParamExtractor.hpp"
#include "ActsExamples/TruthTracking/TrackModifier.hpp"
#include "ActsExamples/TruthTracking/TrackParameterSelector.hpp"
#include "ActsExamples/TruthTracking/TrackParameterSmearing.hpp"
#include "ActsExamples/TruthTracking/TrackTruthMatcher.hpp"
#include "ActsExamples/TruthTracking/TruthSeedingAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/TruthTracking/TruthVertexFinder.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace ActsPython {

void addTruthTracking(py::module& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(TruthTrackFinder, mex, "TruthTrackFinder",
                                inputParticles, inputParticleMeasurementsMap,
                                inputMeasurements, inputSimHits,
                                inputMeasurementSimHitsMap, outputProtoTracks);

  ACTS_PYTHON_DECLARE_ALGORITHM(ParticleTrackParamExtractor, mex,
                                "ParticleTrackParamExtractor", inputParticles,
                                outputTrackParameters);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackParameterSmearing, mex, "TrackParameterSmearing",
      inputTrackParameters, outputTrackParameters, sigmaLoc0, sigmaLoc0PtA,
      sigmaLoc0PtB, sigmaLoc1, sigmaLoc1PtA, sigmaLoc1PtB, sigmaTime, sigmaPhi,
      sigmaTheta, sigmaPtRel, initialSigmas, initialSigmaQoverPt,
      initialSigmaPtRel, initialVarInflation, particleHypothesis,
      randomNumbers);

  {
    using Alg = ParticleSelector;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, IAlgorithm, std::shared_ptr<Alg>>(
                   mex, "ParticleSelector")
                   .def(py::init<const Alg::Config&, Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    {
      auto mc = py::class_<Alg::MeasurementCounter>(alg, "MeasurementCounter")
                    .def(py::init<>())
                    .def("addCounter", &Alg::MeasurementCounter::addCounter);
    }

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(
        c, inputParticles, inputParticleMeasurementsMap, inputMeasurements,
        outputParticles, rhoMin, rhoMax, absZMin, absZMax, timeMin, timeMax,
        phiMin, phiMax, etaMin, etaMax, absEtaMin, absEtaMax, mMin, mMax, ptMin,
        ptMax, hitsMin, hitsMax, measurementsMin, measurementsMax,
        removeCharged, removeNeutral, removeSecondaries, excludeAbsPdgs,
        minPrimaryVertexId, maxPrimaryVertexId, measurementCounter);

    pythonRangeProperty(c, "rho", &Config::rhoMin, &Config::rhoMax);
    pythonRangeProperty(c, "absZ", &Config::absZMin, &Config::absZMax);
    pythonRangeProperty(c, "time", &Config::timeMin, &Config::timeMax);
    pythonRangeProperty(c, "phi", &Config::phiMin, &Config::phiMax);
    pythonRangeProperty(c, "eta", &Config::etaMin, &Config::etaMax);
    pythonRangeProperty(c, "absEta", &Config::absEtaMin, &Config::absEtaMax);
    pythonRangeProperty(c, "m", &Config::mMin, &Config::mMax);
    pythonRangeProperty(c, "pt", &Config::ptMin, &Config::ptMax);
    pythonRangeProperty(c, "measurements", &Config::measurementsMin,
                        &Config::measurementsMax);
    pythonRangeProperty(c, "hits", &Config::hitsMin, &Config::hitsMax);
    pythonRangeProperty(c, "primaryVertexId", &Config::minPrimaryVertexId,
                        &Config::maxPrimaryVertexId);
  }

  {
    using Alg = TrackParameterSelector;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, IAlgorithm, std::shared_ptr<Alg>>(
                   mex, "TrackParameterSelector")
                   .def(py::init<const Alg::Config&, Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, inputTrackParameters, outputTrackParameters, loc0Min,
                       loc0Max, loc1Min, loc1Max, timeMin, timeMax, phiMin,
                       phiMax, etaMin, etaMax, absEtaMin, absEtaMax, ptMin,
                       ptMax);

    pythonRangeProperty(c, "loc0", &Config::loc0Min, &Config::loc0Max);
    pythonRangeProperty(c, "loc1", &Config::loc1Min, &Config::loc1Max);
    pythonRangeProperty(c, "time", &Config::timeMin, &Config::timeMax);
    pythonRangeProperty(c, "phi", &Config::phiMin, &Config::phiMax);
    pythonRangeProperty(c, "eta", &Config::etaMin, &Config::etaMax);
    pythonRangeProperty(c, "absEta", &Config::absEtaMin, &Config::absEtaMax);
    pythonRangeProperty(c, "pt", &Config::ptMin, &Config::ptMax);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(TruthVertexFinder, mex, "TruthVertexFinder",
                                inputTracks, inputTrackParticleMatching,
                                outputProtoVertices, excludeSecondaries,
                                separateSecondaries);

  ACTS_PYTHON_DECLARE_ALGORITHM(TrackModifier, mex, "TrackModifier",
                                inputTracks, outputTracks, dropCovariance,
                                covScale, killTime);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TruthSeedingAlgorithm, mex, "TruthSeedingAlgorithm", inputParticles,
      inputParticleMeasurementsMap, inputSimHits, inputMeasurementSimHitsMap,
      inputSpacePoints, outputParticles, outputSeeds, outputProtoTracks,
      outputParticleHypotheses, particleHypothesis, deltaRMin, deltaRMax,
      absDeltaZMin, absDeltaZMax);

  ACTS_PYTHON_DECLARE_ALGORITHM(HitSelector, mex, "HitSelector", inputHits,
                                inputParticlesSelected, outputHits, minX, maxX,
                                minY, maxY, minZ, maxZ, minR, maxR, minTime,
                                maxTime, minEnergyLoss, maxEnergyLoss,
                                minPrimaryVertexId, maxPrimaryVertexId);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackTruthMatcher, mex, "TrackTruthMatcher", inputTracks, inputParticles,
      inputMeasurementParticlesMap, outputTrackParticleMatching,
      outputParticleTrackMatching, matchingRatio, doubleMatching);
}

}  // namespace ActsPython
