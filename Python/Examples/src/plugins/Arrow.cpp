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
#include "ActsExamples/Io/Arrow/ColliderMLRelease1InputConverter.hpp"
#include "ActsExamples/Io/Parquet/ArrowOutputConverter.hpp"
#include "ActsExamples/Io/Parquet/ParquetReader.hpp"
#include "ActsExamples/Io/Parquet/ParquetWriter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsPython;
using namespace ActsExamples;

PYBIND11_MODULE(ActsExamplesPythonBindingsArrow, m) {
  // ArrowSchemaHandle (Python: `acts.arrow.ArrowSchema`) is registered by
  // the plugin-level binding `ActsPluginsPythonBindingsArrow`. The
  // generated `acts/examples/arrow.py` does `from acts import arrow`
  // before importing this module so that the type is in pybind's
  // registry by the time `expectedSchemas` is bound below.
  ACTS_PYTHON_DECLARE_READER(ParquetReader, m, "ParquetReader", inputDir,
                             collections, expectedSchemas);

  ACTS_PYTHON_DECLARE_WRITER(ParquetWriter, m, "ParquetWriter", outputDir,
                             collections, expectedSchemas, eventsPerShard,
                             eventsPerRowGroup, maxOpenShards);

  py::class_<ArrowOutputConverter, IAlgorithm,
             std::shared_ptr<ArrowOutputConverter>>(m, "ArrowOutputConverter")
      .def_property_readonly("collections", &ArrowOutputConverter::collections);

  {
    auto [alg, c] =
        declareAlgorithm<ArrowParticleOutputConverter, ArrowOutputConverter>(
            m, "ArrowParticleOutputConverter");
    ACTS_PYTHON_STRUCT(c, inputParticles, outputTable);
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
    ACTS_PYTHON_STRUCT(c, inputSimHits, inputParticles, inputClusters,
                       inputSimHitMeasurementsMap, outputTable,
                       detectorResolver);
  }

  {
    auto [alg, c] = declareAlgorithm<ColliderMLRelease1InputConverter,
                                     ActsExamples::IAlgorithm>(
        m, "ColliderMLRelease1InputConverter");
    ACTS_PYTHON_STRUCT(c, inputParticlesTable, inputHitsTable, outputParticles,
                       outputSimHits, outputMeasurements, outputClusters,
                       outputMeasurementSubset, outputMeasSimHitsMap,
                       outputMeasParticlesMap, outputParticleMeasurementsMap,
                       trackingGeometry, geoIdMapPath, geoIdMapSourcePrefix,
                       geoIdMapTargetPrefix, hitBoundsTolerance,
                       keepParticlesWithoutHits);

    alg.def_static(
        "particleSchema", &ColliderMLRelease1InputConverter::particleSchema,
        "Expected schema for the ColliderML Release 1 per-event particle "
        "table.");
    alg.def_static(
        "hitSchema", &ColliderMLRelease1InputConverter::hitSchema,
        "Expected schema for the ColliderML Release 1 per-event tracker-hit "
        "table.");
  }

  m.def("makeVolumeIdDetectorResolver",
        &ArrowSimHitOutputConverter::makeVolumeIdDetectorResolver,
        "volumeToDetector"_a, "defaultValue"_a = static_cast<std::uint8_t>(255),
        R"doc(
Create a C++ detector resolver from a volume-id lookup map.

This avoids per-hit Python callback overhead by baking the mapping into
a C++ lambda once at configuration time.
)doc");
}
