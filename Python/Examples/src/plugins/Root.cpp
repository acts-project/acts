// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/Io/Root/RootAthenaDumpReader.hpp"
#include "ActsExamples/Io/Root/RootAthenaNTupleReader.hpp"
#include "ActsExamples/Io/Root/RootBFieldWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackReader.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialWriter.hpp"
#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"
#include "ActsExamples/Io/Root/RootMuonSpacePointReader.hpp"
#include "ActsExamples/Io/Root/RootMuonSpacePointWriter.hpp"
#include "ActsExamples/Io/Root/RootNuclearInteractionParametersWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"
#include "ActsExamples/Io/Root/RootPropagationSummaryWriter.hpp"
#include "ActsExamples/Io/Root/RootSeedWriter.hpp"
#include "ActsExamples/Io/Root/RootSimHitReader.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/Io/Root/RootSpacepointWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackFinderNTupleWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryReader.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryWriter.hpp"
#include "ActsExamples/Io/Root/RootVertexNTupleWriter.hpp"
#include "ActsExamples/Io/Root/RootVertexReader.hpp"
#include "ActsExamples/Io/Root/RootVertexWriter.hpp"
#include "ActsExamples/Root/MuonVisualization.hpp"
#include "ActsExamples/Root/ScalingCalibrator.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsRoot, root) {
  // Input
  {
    ACTS_PYTHON_DECLARE_READER(RootParticleReader, root, "RootParticleReader",
                               outputParticles, treeName, filePath);

    ACTS_PYTHON_DECLARE_READER(RootVertexReader, root, "RootVertexReader",
                               outputVertices, treeName, filePath);

    ACTS_PYTHON_DECLARE_READER(
        RootMaterialTrackReader, root, "RootMaterialTrackReader",
        outputMaterialTracks, treeName, fileList, readCachedSurfaceInformation);

    ACTS_PYTHON_DECLARE_READER(RootTrackSummaryReader, root,
                               "RootTrackSummaryReader", outputTracks,
                               outputParticles, treeName, filePath);
    ACTS_PYTHON_DECLARE_READER(RootMuonSpacePointReader, root,
                               "RootMuonSpacePointReader", outputSpacePoints,
                               filePath, treeName);

    ACTS_PYTHON_DECLARE_READER(
        RootAthenaNTupleReader, root, "RootAthenaNTupleReader", inputTreeName,
        inputFilePath, outputTrackParameters, outputTruthVtxParameters,
        outputRecoVtxParameters, outputBeamspotConstraint);

    ACTS_PYTHON_DECLARE_READER(
        RootAthenaDumpReader, root, "RootAthenaDumpReader", treename,
        inputfiles, outputMeasurements, outputPixelSpacePoints,
        outputStripSpacePoints, outputSpacePoints, outputClusters,
        outputMeasurementParticlesMap, outputParticleMeasurementsMap,
        outputParticles, onlySpacepoints, onlyPassedParticles,
        skipOverlapSPsPhi, skipOverlapSPsEta, geometryIdMap, trackingGeometry,
        absBoundaryTolerance, onlySpacepoints, noTruth, readCellData);

#ifdef WITH_GEOMODEL_PLUGIN
    ACTS_PYTHON_DECLARE_READER(RootAthenaDumpGeoIdCollector, root,
                               "RootAthenaDumpGeoIdCollector", treename,
                               inputfile, trackingGeometry, geometryIdMap);
#endif

    ACTS_PYTHON_DECLARE_READER(RootSimHitReader, root, "RootSimHitReader",
                               treeName, filePath, outputSimHits);
  }

  // Output
  {
    // Bindings for axis variant used in PlotTool configurations
    {
      py::class_<Acts::Experimental::AxisVariant>(root, "AxisVariant")
          .def_static(
              "regular",
              [](unsigned bins, double low, double high,
                 const std::string& title) {
                return Acts::Experimental::AxisVariant(
                    Acts::Experimental::BoostRegularAxis(bins, low, high,
                                                         title));
              },
              "bins"_a, "low"_a, "high"_a, "title"_a = "")
          .def_static(
              "log",
              [](unsigned bins, double low, double high,
                 const std::string& title) {
                return Acts::Experimental::AxisVariant(
                    Acts::Experimental::BoostLogAxis(bins, low, high, title));
              },
              "bins"_a, "low"_a, "high"_a, "title"_a = "")
          .def_static(
              "variable",
              [](const std::vector<double>& edges, const std::string& title) {
                return Acts::Experimental::AxisVariant(
                    Acts::Experimental::BoostVariableAxis(edges, title));
              },
              "edges"_a, "title"_a = "");

      py::class_<EffPlotTool::Config>(root, "EffPlotToolConfig")
          .def(py::init<>())
          .def_readwrite("varBinning", &EffPlotTool::Config::varBinning)
          .def_readwrite("minTruthPt", &EffPlotTool::Config::minTruthPt);

      py::class_<FakePlotTool::Config>(root, "FakePlotToolConfig")
          .def(py::init<>())
          .def_readwrite("varBinning", &FakePlotTool::Config::varBinning);

      py::class_<DuplicationPlotTool::Config>(root, "DuplicationPlotToolConfig")
          .def(py::init<>())
          .def_readwrite("varBinning",
                         &DuplicationPlotTool::Config::varBinning);

      py::class_<ResPlotTool::Config>(root, "ResPlotToolConfig")
          .def(py::init<>())
          .def_readwrite("varBinning", &ResPlotTool::Config::varBinning);

      py::class_<TrackQualityPlotTool::Config>(root,
                                               "TrackQualityPlotToolConfig")
          .def(py::init<>())
          .def_readwrite("varBinning",
                         &TrackQualityPlotTool::Config::varBinning);

      py::class_<TrackSummaryPlotTool::Config>(root,
                                               "TrackSummaryPlotToolConfig")
          .def(py::init<>())
          .def_readwrite("varBinning",
                         &TrackSummaryPlotTool::Config::varBinning);
    }

    // ROOT WRITERS
    ACTS_PYTHON_DECLARE_WRITER(RootPropagationStepsWriter, root,
                               "RootPropagationStepsWriter", collection,
                               filePath, fileMode);

    ACTS_PYTHON_DECLARE_WRITER(RootPropagationSummaryWriter, root,
                               "RootPropagationSummaryWriter",
                               inputSummaryCollection, filePath, fileMode);

    ACTS_PYTHON_DECLARE_WRITER(RootParticleWriter, root, "RootParticleWriter",
                               inputParticles, filePath, fileMode, treeName,
                               referencePoint, bField, writeHelixParameters);

    ACTS_PYTHON_DECLARE_WRITER(RootVertexWriter, root, "RootVertexWriter",
                               inputVertices, filePath, fileMode, treeName);

    ACTS_PYTHON_DECLARE_WRITER(RootMuonSpacePointWriter, root,
                               "RootMuonSpacePointWriter", inputSpacePoints,
                               filePath, fileMode, treeName, trackingGeometry,
                               writeGlobal);

    ACTS_PYTHON_DECLARE_WRITER(RootTrackFinderNTupleWriter, root,
                               "RootTrackFinderNTupleWriter", inputTracks,
                               inputParticles, inputParticleMeasurementsMap,
                               inputTrackParticleMatching, filePath, fileMode,
                               treeNameTracks, treeNameParticles);

    ACTS_PYTHON_DECLARE_WRITER(RootTrackFitterPerformanceWriter, root,
                               "RootTrackFitterPerformanceWriter", inputTracks,
                               inputParticles, inputTrackParticleMatching,
                               filePath, resPlotToolConfig, effPlotToolConfig,
                               trackSummaryPlotToolConfig);

    ACTS_PYTHON_DECLARE_WRITER(
        RootTrackParameterWriter, root, "RootTrackParameterWriter",
        inputTrackParameters, inputProtoTracks, inputParticles, inputSimHits,
        inputMeasurementParticlesMap, inputMeasurementSimHitsMap, filePath,
        treeName, fileMode);

    ACTS_PYTHON_DECLARE_WRITER(
        RootMaterialTrackWriter, root, "RootMaterialTrackWriter",
        inputMaterialTracks, filePath, fileMode, treeName, recalculateTotals,
        prePostStep, storeSurface, storeVolume, collapseInteractions);

    {
      using Writer = RootBFieldWriter;
      auto w = py::class_<Writer>(root, "RootBFieldWriter")
                   .def_static(
                       "run",
                       [](const Writer::Config& config, Logging::Level level) {
                         Writer::run(config, getDefaultLogger(
                                                 "RootBFieldWriter", level));
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
                   root, "RootMeasurementWriter")
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
              root, "RootMaterialWriter")
              .def(py::init<const Writer::Config&, Logging::Level>(),
                   py::arg("config"), py::arg("level"))
              .def("write",
                   py::overload_cast<const TrackingGeometry&>(&Writer::write));

      auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

      ACTS_PYTHON_STRUCT(c, processSensitives, processApproaches,
                         processRepresenting, processBoundaries, accessorConfig,
                         accessorOptions, filePath, fileMode);
    }

    ACTS_PYTHON_DECLARE_WRITER(RootSeedWriter, root, "RootSeedWriter",
                               inputSeeds, writingMode, filePath, fileMode,
                               treeName);

    ACTS_PYTHON_DECLARE_WRITER(RootSimHitWriter, root, "RootSimHitWriter",
                               inputSimHits, filePath, fileMode, treeName);

    ACTS_PYTHON_DECLARE_WRITER(
        RootSpacepointWriter, root, "RootSpacepointWriter", inputSpacepoints,
        inputMeasurementParticlesMap, filePath, fileMode, treeName);

    ACTS_PYTHON_DECLARE_WRITER(
        RootTrackStatesWriter, root, "RootTrackStatesWriter", inputTracks,
        inputParticles, inputTrackParticleMatching, inputSimHits,
        inputMeasurementSimHitsMap, filePath, treeName, fileMode);

    ACTS_PYTHON_DECLARE_WRITER(
        RootTrackSummaryWriter, root, "RootTrackSummaryWriter", inputTracks,
        inputParticles, inputTrackParticleMatching, filePath, treeName,
        fileMode, writeCovMat, writeGsfSpecific, writeGx2fSpecific);

    ACTS_PYTHON_DECLARE_WRITER(
        RootVertexNTupleWriter, root, "RootVertexNTupleWriter", inputVertices,
        inputTracks, inputTruthVertices, inputParticles, inputSelectedParticles,
        inputTrackParticleMatching, inputVertexTruthMatching, bField, filePath,
        treeName, fileMode, writeTrackInfo);

    ACTS_PYTHON_DECLARE_WRITER(
        RootTrackFinderPerformanceWriter, root,
        "RootTrackFinderPerformanceWriter", inputTracks, inputParticles,
        inputTrackParticleMatching, inputParticleTrackMatching,
        inputParticleMeasurementsMap, filePath, fileMode, effPlotToolConfig,
        fakePlotToolConfig, duplicationPlotToolConfig,
        trackSummaryPlotToolConfig, trackQualityPlotToolConfig,
        subDetectorTrackSummaryVolumes, writeMatchingDetails);

    ACTS_PYTHON_DECLARE_WRITER(RootNuclearInteractionParametersWriter, root,
                               "RootNuclearInteractionParametersWriter",
                               inputSimulationProcesses, filePath, fileMode,
                               interactionProbabilityBins, momentumBins,
                               invariantMassBins, multiplicityMax,
                               writeOptionalHistograms, nSimulatedEvents);
  }

  // Calibration
  {
    root.def(
        "makeScalingCalibrator",
        [](const char* path) -> std::shared_ptr<MeasurementCalibrator> {
          return std::make_shared<ScalingCalibrator>(path);
        },
        py::arg("path"));
  }

  // Muon visualization
  {
    root.def("makeMuonVisualizationFunction", []() {
      return std::function<void(
          const std::string&, const Acts::GeometryContext&,
          const MuonSpacePointBucket&, const SimHitContainer&,
          const SimParticleContainer&, const Acts::TrackingGeometry&,
          const Acts::Logger&)>(visualizeMuonSpacePoints);
    });
  }

  // Muon Hough visualization
  {
    root.def("makeMuonHoughVisualizationFunction", []() {
      return std::function<void(
          const std::string&, const MuonSpacePoint::MuonId&,
          const std::vector<
              Acts::HoughTransformUtils::PeakFinders::IslandsAroundMax<
                  const MuonSpacePoint*>::Maximum>&,
          const Acts::HoughTransformUtils::HoughPlane<const MuonSpacePoint*>&,
          const Acts::HoughTransformUtils::HoughAxisRanges&,
          const MuonSegmentContainer&, const Acts::Logger&)>(
          visualizeMuonHoughMaxima);
    });
  }
}
