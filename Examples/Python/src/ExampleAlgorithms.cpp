// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/TrackFinding/TrackSelector.hpp"
#include "ActsExamples/Fatras/FatrasSimulation.hpp"
#include "ActsExamples/Io/Json/JsonGeometryList.hpp"
#include "ActsExamples/Printers/HitsPrinter.hpp"
#include "ActsExamples/Printers/ParticlesPrinter.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Utilities/TrackSelectorAlgorithm.hpp"

#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addExampleAlgorithms(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  mex.def("readJsonGeometryList", ActsExamples::readJsonGeometryList);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::FatrasSimulation, mex, "FatrasSimulation", inputParticles,
      outputParticlesInitial, outputParticlesFinal, outputSimHits,
      imputParametrisationNuclearInteraction, randomNumbers, trackingGeometry,
      magneticField, pMin, emScattering, emEnergyLossIonisation,
      emEnergyLossRadiation, emPhotonConversion, generateHitsOnSensitive,
      generateHitsOnMaterial, generateHitsOnPassive, averageHitsPerParticle);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::ParticlesPrinter, mex,
                                "ParticlesPrinter", inputParticles);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::HitsPrinter, mex, "HitsPrinter", inputClusters,
      inputMeasurementParticlesMap, inputHitIds, selectIndexStart,
      selectIndexLength, selectVolume, selectLayer, selectModule);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackParametersPrinter, mex,
                                "TrackParametersPrinter", inputTrackParameters);

  {
    using Alg = ActsExamples::TrackSelectorAlgorithm;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, IAlgorithm, std::shared_ptr<Alg>>(
                   mex, "TrackSelectorAlgorithm")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputTracks);
    ACTS_PYTHON_MEMBER(outputTracks);
    ACTS_PYTHON_MEMBER(selectorConfig);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Config = Acts::TrackSelector::Config;
    using CutConfig = Acts::TrackSelector::CutConfig;

    auto tool = py::class_<Acts::TrackSelector>(m, "TrackSelector")
                    .def(py::init<const Config&>(), py::arg("config"));

    {
      auto c = py::class_<CutConfig>(tool, "CutConfig").def(py::init<>());

      patchKwargsConstructor(c);

      ACTS_PYTHON_STRUCT_BEGIN(c, CutConfig);
      ACTS_PYTHON_MEMBER(loc0Min);
      ACTS_PYTHON_MEMBER(loc0Max);
      ACTS_PYTHON_MEMBER(loc1Min);
      ACTS_PYTHON_MEMBER(loc1Max);
      ACTS_PYTHON_MEMBER(timeMin);
      ACTS_PYTHON_MEMBER(timeMax);
      ACTS_PYTHON_MEMBER(phiMin);
      ACTS_PYTHON_MEMBER(phiMax);
      ACTS_PYTHON_MEMBER(etaMin);
      ACTS_PYTHON_MEMBER(etaMax);
      ACTS_PYTHON_MEMBER(absEtaMin);
      ACTS_PYTHON_MEMBER(absEtaMax);
      ACTS_PYTHON_MEMBER(ptMin);
      ACTS_PYTHON_MEMBER(ptMax);
      ACTS_PYTHON_MEMBER(minMeasurements);
      ACTS_PYTHON_STRUCT_END();

      pythonRangeProperty(c, "loc0", &CutConfig::loc0Min, &CutConfig::loc0Max);
      pythonRangeProperty(c, "loc1", &CutConfig::loc1Min, &CutConfig::loc1Max);
      pythonRangeProperty(c, "time", &CutConfig::timeMin, &CutConfig::timeMax);
      pythonRangeProperty(c, "phi", &CutConfig::phiMin, &CutConfig::phiMax);
      pythonRangeProperty(c, "eta", &CutConfig::etaMin, &CutConfig::etaMax);
      pythonRangeProperty(c, "absEta", &CutConfig::absEtaMin,
                          &CutConfig::absEtaMax);
      pythonRangeProperty(c, "pt", &CutConfig::ptMin, &CutConfig::ptMax);
    }

    {
      auto c = py::class_<Config>(tool, "Config")
                   .def(py::init<>())
                   .def(py::init<const CutConfig&>());

      c.def_property_readonly("nEtaBins", &Config::nEtaBins);

      ACTS_PYTHON_STRUCT_BEGIN(c, Config);
      ACTS_PYTHON_MEMBER(cutSets);
      ACTS_PYTHON_MEMBER(absEtaEdges);
      ACTS_PYTHON_STRUCT_END();
    }
  }
}
}  // namespace Acts::Python
