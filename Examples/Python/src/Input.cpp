// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Csv/CsvSpacePointReader.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackReader.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsExamples/Io/Root/RootTrajectorySummaryReader.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace Acts::Python {
void addInput(Context& ctx) {
  auto mex = ctx.get("examples");

  // ROOT READERS
  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootParticleReader, mex,
                             "RootParticleReader", particleCollection,
                             vertexPrimaryCollection, vertexSecondaryCollection,
                             treeName, filePath, orderedEvents);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootMaterialTrackReader, mex,
                             "RootMaterialTrackReader", collection, treeName,
                             fileList, orderedEvents,
                             readCachedSurfaceInformation);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootTrajectorySummaryReader, mex,
                             "RootTrajectorySummaryReader", outputTracks,
                             outputParticles, treeName, filePath,
                             orderedEvents);

  // CSV READERS
  ACTS_PYTHON_DECLARE_READER(ActsExamples::CsvParticleReader, mex,
                             "CsvParticleReader", inputDir, inputStem,
                             outputParticles);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::CsvMeasurementReader, mex,
                             "CsvMeasurementReader", inputDir,
                             outputMeasurements, outputMeasurementSimHitsMap,
                             outputSourceLinks, outputClusters);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::CsvPlanarClusterReader, mex,
                             "CsvPlanarClusterReader", inputDir, outputClusters,
                             outputHitIds, outputMeasurementParticlesMap,
                             outputSimHits, trackingGeometry);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::CsvSimHitReader, mex,
                             "CsvSimHitReader", inputDir, inputStem,
                             outputSimHits);

  ACTS_PYTHON_DECLARE_READER(
      ActsExamples::CsvSpacePointReader, mex, "CsvSpacePointReader", inputDir,
      inputStem, inputCollection, outputSpacePoints, extendCollection);
}
}  // namespace Acts::Python
