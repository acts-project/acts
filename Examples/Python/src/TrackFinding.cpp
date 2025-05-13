// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include "Acts/Seeding/SeedFinderOrthogonalConfig.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SpacePointContainer.hpp"
#include "ActsExamples/TrackFinding/GbtsSeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/HoughTransformSeeder.hpp"
#include "ActsExamples/TrackFinding/MuonHoughSeeder.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedingOrthogonalAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsLookupEstimation.hpp"

#include <cstddef>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addTrackFinding(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::SpacePointMaker, mex, "SpacePointMaker", inputMeasurements,
      outputSpacePoints, trackingGeometry, geometrySelection);

  {
    using Config = Acts::SeedFilterConfig;
    auto c = py::class_<Config>(m, "SeedFilterConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(
        c, deltaInvHelixDiameter, impactWeightFactor, zOriginWeightFactor,
        compatSeedWeight, deltaRMin, maxSeedsPerSpM, compatSeedLimit,
        seedConfirmation, centralSeedConfirmationRange,
        forwardSeedConfirmationRange, useDeltaRorTopRadius, seedWeightIncrement,
        numSeedIncrement, maxSeedsPerSpMConf, maxQualitySeedsPerSpMConf);
    patchKwargsConstructor(c);
  }

  {
    using Config = Acts::SeedFinderConfig<typename Acts::SpacePointContainer<
        ActsExamples::SpacePointContainer<std::vector<const SimSpacePoint*>>,
        Acts::detail::RefHolder>::SpacePointProxyType>;
    auto c = py::class_<Config>(m, "SeedFinderConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(
        c, minPt, cotThetaMax, deltaRMin, deltaRMax, deltaRMinBottomSP,
        deltaRMaxBottomSP, deltaRMinTopSP, deltaRMaxTopSP, impactMax,
        sigmaScattering, maxPtScattering, maxSeedsPerSpM, collisionRegionMin,
        collisionRegionMax, phiMin, phiMax, zMin, zMax, rMax, rMin,
        radLengthPerSeed, zAlign, rAlign, sigmaError, maxBlockSize,
        nTrplPerSpBLimit, nAvgTrplPerSpBLimit, deltaZMax, zBinEdges,
        interactionPointCut, zBinsCustomLooping, useVariableMiddleSPRange,
        deltaRMiddleMinSPRange, deltaRMiddleMaxSPRange, rRangeMiddleSP,
        rMinMiddle, rMaxMiddle, binSizeR, seedConfirmation,
        centralSeedConfirmationRange, forwardSeedConfirmationRange,
        useDetailedDoubleMeasurementInfo);
    patchKwargsConstructor(c);
  }
  {
    using seedOptions = Acts::SeedFinderOptions;
    auto c = py::class_<seedOptions>(m, "SeedFinderOptions").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, beamPos, bFieldInZ);
    patchKwargsConstructor(c);
  }
  {
    using Config =
        Acts::SeedFinderOrthogonalConfig<typename Acts::SpacePointContainer<
            ActsExamples::SpacePointContainer<
                std::vector<const SimSpacePoint*>>,
            Acts::detail::RefHolder>::SpacePointProxyType>;
    auto c =
        py::class_<Config>(m, "SeedFinderOrthogonalConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(
        c, minPt, cotThetaMax, deltaRMinBottomSP, deltaRMaxBottomSP,
        deltaRMinTopSP, deltaRMaxTopSP, impactMax, deltaZMax, sigmaScattering,
        maxPtScattering, maxSeedsPerSpM, collisionRegionMin, collisionRegionMax,
        phiMin, phiMax, zMin, zMax, rMax, rMin, radLengthPerSeed,
        interactionPointCut, deltaPhiMax, highland, maxScatteringAngle2,
        useVariableMiddleSPRange, deltaRMiddleMinSPRange,
        deltaRMiddleMaxSPRange, rRangeMiddleSP, rMinMiddle, rMaxMiddle,
        seedConfirmation, centralSeedConfirmationRange,
        forwardSeedConfirmationRange);
    patchKwargsConstructor(c);
  }

  {
    using Config = Acts::Experimental::SeedFinderGbtsConfig<SimSpacePoint>;
    auto c = py::class_<Config>(m, "SeedFinderGbtsConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, minPt, sigmaScattering, highland, maxScatteringAngle2,
                       ConnectorInputFile, m_phiSliceWidth, m_nMaxPhiSlice,
                       m_useClusterWidth, m_layerGeometry);
    patchKwargsConstructor(c);
  }

  {
    using seedConf = Acts::SeedConfirmationRangeConfig;
    auto c = py::class_<seedConf>(m, "SeedConfirmationRangeConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, zMinSeedConf, zMaxSeedConf, rMaxSeedConf,
                       nTopForLargeR, nTopForSmallR, seedConfMinBottomRadius,
                       seedConfMaxZOrigin, minImpactSeedConf);
    patchKwargsConstructor(c);
  }

  {
    using Config = Acts::CylindricalSpacePointGridConfig;
    auto c = py::class_<Config>(m, "SpacePointGridConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, minPt, rMax, zMax, zMin, phiMin, phiMax, deltaRMax,
                       cotThetaMax, phiBinDeflectionCoverage, maxPhiBins,
                       impactMax, zBinEdges);
    patchKwargsConstructor(c);
  }
  {
    using Options = Acts::CylindricalSpacePointGridOptions;
    auto c = py::class_<Options>(m, "SpacePointGridOptions").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, bFieldInZ);
    patchKwargsConstructor(c);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::SeedingAlgorithm, mex, "SeedingAlgorithm", inputSpacePoints,
      outputSeeds, seedFilterConfig, seedFinderConfig, seedFinderOptions,
      gridConfig, gridOptions, allowSeparateRMax, zBinNeighborsTop,
      zBinNeighborsBottom, numPhiNeighbors, useExtraCuts);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::SeedingOrthogonalAlgorithm, mex,
                                "SeedingOrthogonalAlgorithm", inputSpacePoints,
                                outputSeeds, seedFilterConfig, seedFinderConfig,
                                seedFinderOptions);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::GbtsSeedingAlgorithm, mex, "GbtsSeedingAlgorithm",
      inputSpacePoints, outputSeeds, seedFinderConfig, seedFinderOptions,
      layerMappingFile, geometrySelection, trackingGeometry, ActsGbtsMap,
      fill_module_csv, inputClusters);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::HoughTransformSeeder, mex, "HoughTransformSeeder",
      inputSpacePoints, outputProtoTracks, trackingGeometry, geometrySelection,
      inputMeasurements, subRegions, nLayers, xMin, xMax, yMin, yMax,
      houghHistSize_x, houghHistSize_y, hitExtend_x, threshold,
      localMaxWindowSize, kA);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::MuonHoughSeeder, mex,
                                "MuonHoughSeeder", inTruthSegments,
                                inSpacePoints, outHoughMax);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TrackParamsEstimationAlgorithm, mex,
      "TrackParamsEstimationAlgorithm", inputSeeds, inputProtoTracks,
      outputTrackParameters, outputSeeds, outputProtoTracks, trackingGeometry,
      magneticField, bFieldMin, initialSigmas, initialSigmaQoverPt,
      initialSigmaPtRel, initialVarInflation, noTimeVarInflation,
      particleHypothesis);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackParamsLookupEstimation, mex,
                                "TrackParamsLookupEstimation", refLayers, bins,
                                inputHits, inputParticles,
                                trackLookupGridWriters);

  {
    using Alg = ActsExamples::TrackFindingAlgorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, ActsExamples::IAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackFindingAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config)
            .def_static("makeTrackFinderFunction",
                        [](std::shared_ptr<const Acts::TrackingGeometry>
                               trackingGeometry,
                           std::shared_ptr<const Acts::MagneticFieldProvider>
                               magneticField,
                           Logging::Level level) {
                          return Alg::makeTrackFinderFunction(
                              std::move(trackingGeometry),
                              std::move(magneticField),
                              *Acts::getDefaultLogger("TrackFinding", level));
                        });

    py::class_<Alg::TrackFinderFunction,
               std::shared_ptr<Alg::TrackFinderFunction>>(
        alg, "TrackFinderFunction");

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, inputMeasurements, inputInitialTrackParameters,
                       inputSeeds, outputTracks, trackingGeometry,
                       magneticField, findTracks, measurementSelectorCfg,
                       trackSelectorCfg, maxSteps, twoWay, reverseSearch,
                       seedDeduplication, stayOnSeed, pixelVolumeIds,
                       stripVolumeIds, maxPixelHoles, maxStripHoles, trimTracks,
                       constrainToVolumeIds, endOfWorldVolumeIds);
  }

  {
    auto constructor =
        [](const std::vector<std::pair<
               GeometryIdentifier,
               std::tuple<std::vector<double>, std::vector<double>,
                          std::vector<double>, std::vector<std::size_t>>>>&
               input) {
          std::vector<std::pair<GeometryIdentifier, MeasurementSelectorCuts>>
              converted;
          converted.reserve(input.size());
          for (const auto& [id, cuts] : input) {
            const auto& [bins, chi2Measurement, chi2Outlier, num] = cuts;
            converted.emplace_back(
                id, MeasurementSelectorCuts{bins, chi2Measurement, num,
                                            chi2Outlier});
          }
          return std::make_unique<MeasurementSelector::Config>(converted);
        };

    py::class_<MeasurementSelectorCuts>(m, "MeasurementSelectorCuts")
        .def(py::init<>())
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<std::size_t>, std::vector<double>>())
        .def_readwrite("etaBins", &MeasurementSelectorCuts::etaBins)
        .def_readwrite("chi2CutOffMeasurement",
                       &MeasurementSelectorCuts::chi2CutOff)
        .def_readwrite("chi2CutOffOutlier",
                       &MeasurementSelectorCuts::chi2CutOffOutlier)
        .def_readwrite("numMeasurementsCutOff",
                       &MeasurementSelectorCuts::numMeasurementsCutOff);

    auto ms = py::class_<MeasurementSelector>(m, "MeasurementSelector");
    auto c =
        py::class_<MeasurementSelector::Config>(ms, "Config")
            .def(py::init<std::vector<
                     std::pair<GeometryIdentifier, MeasurementSelectorCuts>>>())
            .def(py::init(constructor));
  }
}

}  // namespace Acts::Python
