// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootAthenaDumpReader.hpp"
#include "ActsExamples/Io/Root/RootAthenaNTupleReader.hpp"
#include "ActsExamples/Io/Root/RootMaterialDecorator.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackReader.hpp"
#include "ActsExamples/Io/Root/RootMuonSpacePointReader.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsExamples/Io/Root/RootSimHitReader.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryReader.hpp"
#include "ActsExamples/Io/Root/RootVertexReader.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace ActsPython {

void addRootInput(Context& ctx) {
  auto mex = ctx.get("examples");

  ACTS_PYTHON_DECLARE_READER(RootParticleReader, mex, "RootParticleReader",
                             outputParticles, treeName, filePath);

  ACTS_PYTHON_DECLARE_READER(RootVertexReader, mex, "RootVertexReader",
                             outputVertices, treeName, filePath);

  ACTS_PYTHON_DECLARE_READER(RootMaterialTrackReader, mex,
                             "RootMaterialTrackReader", outputMaterialTracks,
                             treeName, fileList, readCachedSurfaceInformation);

  ACTS_PYTHON_DECLARE_READER(RootTrackSummaryReader, mex,
                             "RootTrackSummaryReader", outputTracks,
                             outputParticles, treeName, filePath);
  ACTS_PYTHON_DECLARE_READER(RootMuonSpacePointReader, mex,
                             "RootMuonSpacePointReader", outputSpacePoints,
                             filePath, treeName);

  ACTS_PYTHON_DECLARE_READER(
      RootAthenaNTupleReader, mex, "RootAthenaNTupleReader", inputTreeName,
      inputFilePath, outputTrackParameters, outputTruthVtxParameters,
      outputRecoVtxParameters, outputBeamspotConstraint);

  ACTS_PYTHON_DECLARE_READER(
      RootAthenaDumpReader, mex, "RootAthenaDumpReader", treename, inputfiles,
      outputMeasurements, outputPixelSpacePoints, outputStripSpacePoints,
      outputSpacePoints, outputClusters, outputMeasurementParticlesMap,
      outputParticleMeasurementsMap, outputParticles, onlySpacepoints,
      onlyPassedParticles, skipOverlapSPsPhi, skipOverlapSPsEta, geometryIdMap,
      trackingGeometry, absBoundaryTolerance, onlySpacepoints, noTruth,
      readCellData);

#ifdef WITH_GEOMODEL_PLUGIN
  ACTS_PYTHON_DECLARE_READER(RootAthenaDumpGeoIdCollector, mex,
                             "RootAthenaDumpGeoIdCollector", treename,
                             inputfile, trackingGeometry, geometryIdMap);
#endif

  ACTS_PYTHON_DECLARE_READER(RootSimHitReader, mex, "RootSimHitReader",
                             treeName, filePath, outputSimHits);

  {
    auto rmd =
        py::class_<RootMaterialDecorator, IMaterialDecorator,
                   std::shared_ptr<RootMaterialDecorator>>(
            mex, "RootMaterialDecorator")
            .def(py::init<RootMaterialDecorator::Config, Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    using Config = RootMaterialDecorator::Config;
    auto c = py::class_<Config>(rmd, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, accessorConfig, accessorOptions, fileName);
  }
}

}  // namespace ActsPython
