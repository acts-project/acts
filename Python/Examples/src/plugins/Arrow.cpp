// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Io/Arrow/ArrowParticleOutputConverter.hpp"
#include "ActsExamples/Io/Arrow/ArrowSimHitOutputConverter.hpp"
#include "ActsExamples/Io/Arrow/ArrowTrackOutputConverter.hpp"
#include "ActsExamples/Io/Parquet/ArrowOutputConverter.hpp"
#include "ActsExamples/Io/Parquet/ParquetReader.hpp"
#include "ActsExamples/Io/Parquet/ParquetWriter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsPython;
using namespace ActsExamples;

PYBIND11_MODULE(ActsExamplesPythonBindingsArrow, m) {
  ACTS_PYTHON_DECLARE_READER(ParquetReader, m, "ParquetReader", inputDir,
                             collections);

  ACTS_PYTHON_DECLARE_WRITER(ParquetWriter, m, "ParquetWriter", outputDir,
                             collections, eventsPerRowGroup);

  py::class_<ArrowOutputConverter, IAlgorithm,
             std::shared_ptr<ArrowOutputConverter>>(m, "ArrowOutputConverter")
      .def_property_readonly("collections", &ArrowOutputConverter::collections);

  {
    auto [alg, c] =
        declareAlgorithm<ArrowParticleOutputConverter, ArrowOutputConverter>(
            m, "ArrowParticleOutputConverter");
    ACTS_PYTHON_STRUCT(c, inputParticles, outputTable, referencePoint, bField,
                       writeHelixParameters);
  }

  {
    auto [alg, c] =
        declareAlgorithm<ArrowTrackOutputConverter, ArrowOutputConverter>(
            m, "ArrowTrackOutputConverter");
    ACTS_PYTHON_STRUCT(c, inputTracks, inputTrackParticleMatching,
                       inputParticles, inputMeasurementSimHitsMap, outputTable,
                       writeTime);
  }

  {
    auto [alg, c] =
        declareAlgorithm<ArrowSimHitOutputConverter, ArrowOutputConverter>(
            m, "ArrowSimHitOutputConverter");
    ACTS_PYTHON_STRUCT(c, inputSimHits, inputParticles, inputMeasurements,
                       inputSimHitMeasurementsMap, outputTable,
                       trackingGeometry, subsystemByVolume);
  }
}
