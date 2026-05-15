// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/Framework/BufferedReader.hpp"
#include "ActsExamples/Io/Csv/CsvGnnGraphReader.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"
#include "ActsExamples/Io/Csv/CsvMuonSegmentReader.hpp"
#include "ActsExamples/Io/Csv/CsvMuonSpacePointReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Csv/CsvSpacePointReader.hpp"
#include "ActsExamples/Io/Csv/CsvTrackParameterReader.hpp"
#include "ActsExamples/TrackFinding/ITrackParamsLookupReader.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace ActsPython {

void addInput(py::module& mex) {
  // Buffered reader
  ACTS_PYTHON_DECLARE_READER(BufferedReader, mex, "BufferedReader",
                             upstreamReader, selectionSeed, bufferSize);

  ACTS_PYTHON_DECLARE_READER(CsvParticleReader, mex, "CsvParticleReader",
                             inputDir, inputStem, outputParticles);

  ACTS_PYTHON_DECLARE_READER(CsvMeasurementReader, mex, "CsvMeasurementReader",
                             inputDir, outputMeasurements,
                             outputMeasurementSimHitsMap, outputClusters,
                             outputMeasurementParticlesMap, inputSimHits);

  ACTS_PYTHON_DECLARE_READER(CsvSimHitReader, mex, "CsvSimHitReader", inputDir,
                             inputStem, outputSimHits);
  ACTS_PYTHON_DECLARE_READER(CsvMuonSegmentReader, mex, "CsvMuonSegmentReader",
                             inputDir, inputStem, outputSegments);
  ACTS_PYTHON_DECLARE_READER(CsvMuonSpacePointReader, mex,
                             "CsvMuonSpacePointReader", inputDir, inputStem,
                             outputSpacePoints);

  ACTS_PYTHON_DECLARE_READER(CsvSpacePointReader, mex, "CsvSpacePointReader",
                             inputDir, inputStem, inputCollection,
                             outputSpacePoints, extendCollection);

  ACTS_PYTHON_DECLARE_READER(CsvTrackParameterReader, mex,
                             "CsvTrackParameterReader", inputDir, inputStem,
                             outputTrackParameters, beamspot);

  ACTS_PYTHON_DECLARE_READER(CsvGnnGraphReader, mex, "CsvGnnGraphReader",
                             inputDir, inputStem, outputGraph);

  py::class_<ITrackParamsLookupReader,
             std::shared_ptr<ITrackParamsLookupReader>>(
      mex, "ITrackParamsLookupReader");
}

}  // namespace ActsPython
