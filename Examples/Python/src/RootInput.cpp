// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/Root/RootAthenaDumpReader.hpp"
#include "ActsExamples/Io/Root/RootAthenaNTupleReader.hpp"
#include "ActsExamples/Io/Root/RootMaterialDecorator.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackReader.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsExamples/Io/Root/RootSimHitReader.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryReader.hpp"
#include "ActsExamples/Io/Root/RootVertexReader.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace Acts::Python {

void addRootInput(Context& ctx) {
  auto mex = ctx.get("examples");

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootParticleReader, mex,
                             "RootParticleReader", outputParticles, treeName,
                             filePath);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootVertexReader, mex,
                             "RootVertexReader", outputVertices, treeName,
                             filePath);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootMaterialTrackReader, mex,
                             "RootMaterialTrackReader", outputMaterialTracks,
                             treeName, fileList, readCachedSurfaceInformation);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootTrackSummaryReader, mex,
                             "RootTrackSummaryReader", outputTracks,
                             outputParticles, treeName, filePath);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootAthenaNTupleReader, mex,
                             "RootAthenaNTupleReader", inputTreeName,
                             inputFilePath, outputTrackParameters,
                             outputTruthVtxParameters, outputRecoVtxParameters,
                             outputBeamspotConstraint);

  ACTS_PYTHON_DECLARE_READER(
      ActsExamples::RootAthenaDumpReader, mex, "RootAthenaDumpReader", treename,
      inputfile, outputMeasurements, outputPixelSpacePoints,
      outputStripSpacePoints, outputSpacePoints, outputClusters,
      outputMeasurementParticlesMap, outputParticles, onlyPassedParticles,
      skipOverlapSPsPhi, skipOverlapSPsEta, geometryIdMap, trackingGeometry,
      absBoundaryTolerance);

#ifdef WITH_GEOMODEL_PLUGIN
  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootAthenaDumpGeoIdCollector, mex,
                             "RootAthenaDumpGeoIdCollector", treename,
                             inputfile, trackingGeometry, geometryIdMap);
#endif

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootSimHitReader, mex,
                             "RootSimHitReader", treeName, filePath,
                             outputSimHits);
}

{
  auto rmd =
      py::class_<RootMaterialDecorator, Acts::IMaterialDecorator,
                 std::shared_ptr<RootMaterialDecorator>>(
          mex, "RootMaterialDecorator")
          .def(py::init<RootMaterialDecorator::Config, Acts::Logging::Level>(),
               py::arg("config"), py::arg("level"));

  using Config = RootMaterialDecorator::Config;
  auto c = py::class_<Config>(rmd, "Config").def(py::init<>());

  ACTS_PYTHON_STRUCT_BEGIN(c, Config);
  ACTS_PYTHON_MEMBER(voltag);
  ACTS_PYTHON_MEMBER(boutag);
  ACTS_PYTHON_MEMBER(laytag);
  ACTS_PYTHON_MEMBER(apptag);
  ACTS_PYTHON_MEMBER(sentag);
  ACTS_PYTHON_MEMBER(ntag);
  ACTS_PYTHON_MEMBER(vtag);
  ACTS_PYTHON_MEMBER(otag);
  ACTS_PYTHON_MEMBER(mintag);
  ACTS_PYTHON_MEMBER(maxtag);
  ACTS_PYTHON_MEMBER(ttag);
  ACTS_PYTHON_MEMBER(x0tag);
  ACTS_PYTHON_MEMBER(l0tag);
  ACTS_PYTHON_MEMBER(atag);
  ACTS_PYTHON_MEMBER(ztag);
  ACTS_PYTHON_MEMBER(rhotag);
  ACTS_PYTHON_MEMBER(fileName);
  ACTS_PYTHON_STRUCT_END();
}

}  // namespace Acts::Python
