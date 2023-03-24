// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementReader.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementWriter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepMultiTrajectoryWriter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepParticleReader.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepParticleWriter.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitReader.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitWriter.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace Acts::Python {

void addEDM4hep(Context& ctx) {
  auto mex = ctx.get("examples");
  auto edm4hep = mex.def_submodule("_edm4hep");

  ACTS_PYTHON_DECLARE_READER(ActsExamples::EDM4hepSimHitReader, edm4hep,
                             "EDM4hepSimHitReader", inputPath, inputParticles,
                             outputSimHits, dd4hepDetector);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::EDM4hepSimHitWriter, edm4hep,
                             "EDM4hepSimHitWriter", inputSimHits,
                             inputParticles, outputPath, outputParticles,
                             outputSimTrackerHits);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::EDM4hepMeasurementReader, edm4hep,
                             "EDM4hepMeasurementReader", inputPath,
                             outputMeasurements, outputMeasurementSimHitsMap,
                             outputSourceLinks, outputClusters);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::EDM4hepMeasurementWriter, edm4hep,
                             "EDM4hepMeasurementWriter", inputMeasurements,
                             inputClusters, outputPath);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::EDM4hepParticleReader, edm4hep,
                             "EDM4hepParticleReader", inputPath, inputParticles,
                             outputParticles);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::EDM4hepParticleWriter, edm4hep,
                             "EDM4hepParticleWriter", inputParticles,
                             outputPath, outputParticles);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::EDM4hepMultiTrajectoryWriter,
                             edm4hep, "EDM4hepMultiTrajectoryWriter",
                             inputTrajectories, inputMeasurementParticlesMap,
                             outputPath);
}

}  // namespace Acts::Python
