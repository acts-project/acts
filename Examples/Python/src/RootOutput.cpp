// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootBFieldWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialWriter.hpp"
#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"
#include "ActsExamples/Io/Root/RootMuonSpacePointWriter.hpp"
#include "ActsExamples/Io/Root/RootNuclearInteractionParametersWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"
#include "ActsExamples/Io/Root/RootPropagationSummaryWriter.hpp"
#include "ActsExamples/Io/Root/RootSeedWriter.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/Io/Root/RootSpacepointWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryWriter.hpp"
#include "ActsExamples/Io/Root/RootVertexWriter.hpp"
#include "ActsExamples/Io/Root/TrackFinderNTupleWriter.hpp"
#include "ActsExamples/Io/Root/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/TrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/VertexNTupleWriter.hpp"
#include "ActsPlugins/Root/RootMaterialMapIo.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace Acts {
class TrackingGeometry;
namespace detail {
struct Step;
}  // namespace detail
}  // namespace Acts
namespace ActsExamples {
class IWriter;
struct AlgorithmContext;
}  // namespace ActsExamples

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace ActsPython {

void addRootOutput(Context& ctx) {
  auto& mex = ctx.get("examples");

  // Bindings for the binning in e.g., TrackFinderPerformanceWriter
  {
    py::class_<PlotHelpers::Binning>(mex, "Binning")
        .def(py::init<std::string, int, double, double>(), "title"_a, "bins"_a,
             "bMin"_a, "bMax"_a)
        .def(py::init<std::string, std::vector<double>>(), "title"_a, "bins"_a);

    py::class_<EffPlotTool::Config>(mex, "EffPlotToolConfig")
        .def(py::init<std::map<std::string, PlotHelpers::Binning>>(),
             "varBinning"_a);

    py::class_<FakePlotTool::Config>(mex, "FakePlotToolConfig")
        .def(py::init<std::map<std::string, PlotHelpers::Binning>>(),
             "varBinning"_a);

    py::class_<DuplicationPlotTool::Config>(mex, "DuplicationPlotToolConfig")
        .def(py::init<std::map<std::string, PlotHelpers::Binning>>(),
             "varBinning"_a);
  }

  // ROOT WRITERS
  ACTS_PYTHON_DECLARE_WRITER(RootPropagationStepsWriter, mex,
                             "RootPropagationStepsWriter", collection, filePath,
                             fileMode);

  ACTS_PYTHON_DECLARE_WRITER(RootPropagationSummaryWriter, mex,
                             "RootPropagationSummaryWriter",
                             inputSummaryCollection, filePath, fileMode);

  ACTS_PYTHON_DECLARE_WRITER(RootParticleWriter, mex, "RootParticleWriter",
                             inputParticles, filePath, fileMode, treeName,
                             referencePoint, bField, writeHelixParameters);

  ACTS_PYTHON_DECLARE_WRITER(RootVertexWriter, mex, "RootVertexWriter",
                             inputVertices, filePath, fileMode, treeName);

  ACTS_PYTHON_DECLARE_WRITER(RootMuonSpacePointWriter, mex,
                             "RootMuonSpacePointWriter", inputSpacePoints,
                             filePath, fileMode, treeName, trackingGeometry,
                             writeGlobal);

  ACTS_PYTHON_DECLARE_WRITER(
      TrackFinderNTupleWriter, mex, "TrackFinderNTupleWriter", inputTracks,
      inputParticles, inputParticleMeasurementsMap, inputTrackParticleMatching,
      filePath, fileMode, treeNameTracks, treeNameParticles);

  ACTS_PYTHON_DECLARE_WRITER(
      TrackFitterPerformanceWriter, mex, "TrackFitterPerformanceWriter",
      inputTracks, inputParticles, inputTrackParticleMatching, filePath,
      resPlotToolConfig, effPlotToolConfig, trackSummaryPlotToolConfig);

  ACTS_PYTHON_DECLARE_WRITER(
      RootTrackParameterWriter, mex, "RootTrackParameterWriter",
      inputTrackParameters, inputProtoTracks, inputParticles, inputSimHits,
      inputMeasurementParticlesMap, inputMeasurementSimHitsMap, filePath,
      treeName, fileMode);

  ACTS_PYTHON_DECLARE_WRITER(
      RootMaterialTrackWriter, mex, "RootMaterialTrackWriter",
      inputMaterialTracks, filePath, fileMode, treeName, recalculateTotals,
      prePostStep, storeSurface, storeVolume, collapseInteractions);

  {
    using Writer = RootBFieldWriter;
    auto w = py::class_<Writer>(mex, "RootBFieldWriter")
                 .def_static(
                     "run",
                     [](const Writer::Config& config, Logging::Level level) {
                       Writer::run(config,
                                   getDefaultLogger("RootBFieldWriter", level));
                     },
                     py::arg("config"), py::arg("level"));

    py::enum_<Writer::GridType>(w, "GridType")
        .value("rz", Writer::GridType::rz)
        .value("xyz", Writer::GridType::xyz);

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, treeName, fileName, fileMode, bField, gridType,
                       rBounds, zBounds, rBins, zBins, phiBins);
  }

  {
    using Writer = RootMeasurementWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "RootMeasurementWriter")
                 .def(py::init<const Writer::Config&, Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, inputMeasurements, inputClusters, inputSimHits,
                       inputMeasurementSimHitsMap, filePath, fileMode,
                       surfaceByIdentifier);
  }

  {
    using Writer = RootMaterialWriter;
    auto w =
        py::class_<Writer, IMaterialWriter, std::shared_ptr<Writer>>(
            mex, "RootMaterialWriter")
            .def(py::init<const Writer::Config&, Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write",
                 py::overload_cast<const TrackingGeometry&>(&Writer::write));

    auto ac =
        py::class_<ActsPlugins::RootMaterialMapIo::Config>(w, "AccessorConfig")
            .def(py::init<>());

    ACTS_PYTHON_STRUCT(ac, volumePrefix, portalPrefix, layerPrefix,
                       passivePrefix, sensitivePrefix, nBinsHistName,
                       axisDirHistName, axisBoundaryTypeHistName, indexHistName,
                       minRangeHistName, maxRangeHistName, thicknessHistName,
                       x0HistName, l0HistName, aHistName, zHistName,
                       rhoHistName);

    auto ao = py::class_<ActsPlugins::RootMaterialMapIo::Options>(
                  w, "AccessorOptions")
                  .def(py::init<>());

    ACTS_PYTHON_STRUCT(ao, homogeneousMaterialTreeName, indexedMaterialTreeName,
                       folderSurfaceNameBase, folderVolumeNameBase,
                       indexedMaterial);

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, processSensitives, processApproaches,
                       processRepresenting, processBoundaries, accessorConfig,
                       accessorOptions, filePath, fileMode);
  }

  ACTS_PYTHON_DECLARE_WRITER(RootSeedWriter, mex, "RootSeedWriter", inputSeeds,
                             writingMode, filePath, fileMode, treeName);

  ACTS_PYTHON_DECLARE_WRITER(RootSimHitWriter, mex, "RootSimHitWriter",
                             inputSimHits, filePath, fileMode, treeName);

  ACTS_PYTHON_DECLARE_WRITER(RootSpacepointWriter, mex, "RootSpacepointWriter",
                             inputSpacepoints, inputMeasurementParticlesMap,
                             filePath, fileMode, treeName);

  ACTS_PYTHON_DECLARE_WRITER(
      RootTrackStatesWriter, mex, "RootTrackStatesWriter", inputTracks,
      inputParticles, inputTrackParticleMatching, inputSimHits,
      inputMeasurementSimHitsMap, filePath, treeName, fileMode);

  ACTS_PYTHON_DECLARE_WRITER(
      RootTrackSummaryWriter, mex, "RootTrackSummaryWriter", inputTracks,
      inputParticles, inputTrackParticleMatching, filePath, treeName, fileMode,
      writeCovMat, writeGsfSpecific, writeGx2fSpecific);

  ACTS_PYTHON_DECLARE_WRITER(
      VertexNTupleWriter, mex, "VertexNTupleWriter", inputVertices, inputTracks,
      inputTruthVertices, inputParticles, inputSelectedParticles,
      inputTrackParticleMatching, bField, filePath, treeName, fileMode,
      vertexMatchThreshold, trackMatchThreshold, writeTrackInfo);

  ACTS_PYTHON_DECLARE_WRITER(
      TrackFinderPerformanceWriter, mex, "TrackFinderPerformanceWriter",
      inputTracks, inputParticles, inputTrackParticleMatching,
      inputParticleTrackMatching, inputParticleMeasurementsMap, filePath,
      fileMode, effPlotToolConfig, fakePlotToolConfig,
      duplicationPlotToolConfig, trackSummaryPlotToolConfig,
      subDetectorTrackSummaryVolumes, writeMatchingDetails);

  ACTS_PYTHON_DECLARE_WRITER(RootNuclearInteractionParametersWriter, mex,
                             "RootNuclearInteractionParametersWriter",
                             inputSimulationProcesses, filePath, fileMode,
                             interactionProbabilityBins, momentumBins,
                             invariantMassBins, multiplicityMax,
                             writeOptionalHistograms, nSimulatedEvents);
}

}  // namespace ActsPython
