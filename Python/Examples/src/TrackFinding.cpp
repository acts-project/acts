// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include "Acts/Seeding/SeedFinderOrthogonalConfig.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SpacePointContainer.hpp"
#include "ActsExamples/TrackFinding/AdaptiveHoughTransformSeeder.hpp"
#include "ActsExamples/TrackFinding/GbtsSeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/GridTripletSeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/HoughTransformSeeder.hpp"
#include "ActsExamples/TrackFinding/MuonHoughSeeder.hpp"
#include "ActsExamples/TrackFinding/OrthogonalTripletSeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedingOrthogonalAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsLookupEstimation.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addTrackFinding(py::module& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(SpacePointMaker, mex, "SpacePointMaker",
                                inputMeasurements, outputSpacePoints,
                                trackingGeometry, geometrySelection,
                                stripGeometrySelection);

  {
    using Config = Acts::SeedFinderConfig<typename Acts::SpacePointContainer<
        SpacePointContainer<std::vector<const SimSpacePoint*>>,
        Acts::detail::RefHolder>::SpacePointProxyType>;
    auto c = py::class_<Config>(mex, "SeedFinderConfig").def(py::init<>());
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
    using Config =
        Acts::SeedFinderOrthogonalConfig<typename Acts::SpacePointContainer<
            SpacePointContainer<std::vector<const SimSpacePoint*>>,
            Acts::detail::RefHolder>::SpacePointProxyType>;
    auto c =
        py::class_<Config>(mex, "SeedFinderOrthogonalConfig").def(py::init<>());
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
    using Config = Acts::Experimental::SeedFinderGbtsConfig;
    auto c = py::class_<Config>(mex, "SeedFinderGbtsConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, minPt, connectorInputFile, phiSliceWidth,
                       nMaxPhiSlice, lutInputFile);
    patchKwargsConstructor(c);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      SeedingAlgorithm, mex, "SeedingAlgorithm", inputSpacePoints, outputSeeds,
      seedFilterConfig, seedFinderConfig, seedFinderOptions, gridConfig,
      gridOptions, allowSeparateRMax, zBinNeighborsTop, zBinNeighborsBottom,
      numPhiNeighbors, useExtraCuts);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      GridTripletSeedingAlgorithm, mex, "GridTripletSeedingAlgorithm",
      inputSpacePoints, outputSeeds, bFieldInZ, minPt, cotThetaMax, impactMax,
      deltaRMin, deltaRMax, deltaRMinTop, deltaRMaxTop, deltaRMinBottom,
      deltaRMaxBottom, rMin, rMax, zMin, zMax, phiMin, phiMax,
      phiBinDeflectionCoverage, maxPhiBins, zBinNeighborsTop,
      zBinNeighborsBottom, numPhiNeighbors, zBinEdges, zBinsCustomLooping,
      rMinMiddle, rMaxMiddle, useVariableMiddleSPRange, rRangeMiddleSP,
      deltaRMiddleMinSPRange, deltaRMiddleMaxSPRange, deltaZMin, deltaZMax,
      interactionPointCut, collisionRegionMin, collisionRegionMax,
      helixCutTolerance, sigmaScattering, radLengthPerSeed, toleranceParam,
      deltaInvHelixDiameter, compatSeedWeight, impactWeightFactor,
      zOriginWeightFactor, maxSeedsPerSpM, compatSeedLimit, seedWeightIncrement,
      numSeedIncrement, seedConfirmation, centralSeedConfirmationRange,
      forwardSeedConfirmationRange, maxSeedsPerSpMConf,
      maxQualitySeedsPerSpMConf, useDeltaRinsteadOfTopRadius, useExtraCuts);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      OrthogonalTripletSeedingAlgorithm, mex,
      "OrthogonalTripletSeedingAlgorithm", inputSpacePoints, outputSeeds,
      bFieldInZ, minPt, cotThetaMax, impactMax, deltaRMin, deltaRMax,
      deltaRMinTop, deltaRMaxTop, deltaRMinBottom, deltaRMaxBottom, rMin, rMax,
      zMin, zMax, phiMin, phiMax, rMinMiddle, rMaxMiddle,
      useVariableMiddleSPRange, rRangeMiddleSP, deltaRMiddleMinSPRange,
      deltaRMiddleMaxSPRange, deltaZMin, deltaZMax, interactionPointCut,
      collisionRegionMin, collisionRegionMax, helixCutTolerance,
      sigmaScattering, radLengthPerSeed, toleranceParam, deltaInvHelixDiameter,
      compatSeedWeight, impactWeightFactor, zOriginWeightFactor, maxSeedsPerSpM,
      compatSeedLimit, seedWeightIncrement, numSeedIncrement, seedConfirmation,
      centralSeedConfirmationRange, forwardSeedConfirmationRange,
      maxSeedsPerSpMConf, maxQualitySeedsPerSpMConf,
      useDeltaRinsteadOfTopRadius, useExtraCuts);

  ACTS_PYTHON_DECLARE_ALGORITHM(SeedingOrthogonalAlgorithm, mex,
                                "SeedingOrthogonalAlgorithm", inputSpacePoints,
                                outputSeeds, seedFilterConfig, seedFinderConfig,
                                seedFinderOptions);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::GbtsSeedingAlgorithm, mex, "GbtsSeedingAlgorithm",
      inputSpacePoints, outputSeeds, seedFinderConfig, seedFinderOptions,
      layerMappingFile, trackingGeometry, actsGbtsMap, fill_module_csv,
      inputClusters);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      HoughTransformSeeder, mex, "HoughTransformSeeder", inputSpacePoints,
      outputProtoTracks, trackingGeometry, geometrySelection, inputMeasurements,
      subRegions, nLayers, xMin, xMax, yMin, yMax, houghHistSize_x,
      houghHistSize_y, hitExtend_x, threshold, localMaxWindowSize, kA);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      AdaptiveHoughTransformSeeder, mex, "AdaptiveHoughTransformSeeder",
      inputSpacePoints, outputSeeds, outputProtoTracks, trackingGeometry,
      qOverPtMin, qOverPtMinBinSize, phiMinBinSize, threshold, noiseThreshold,
      deduplicate, inverseA, doSecondPhase, zRange, cotThetaRange,
      cotThetaMinBinSize, zMinBinSize);

  ACTS_PYTHON_DECLARE_ALGORITHM(MuonHoughSeeder, mex, "MuonHoughSeeder",
                                inTruthSegments, inSpacePoints, outHoughMax,
                                nBinsTanTheta, nBinsY0, nBinsTanPhi, nBinsX0);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackParamsEstimationAlgorithm, mex, "TrackParamsEstimationAlgorithm",
      inputSeeds, inputProtoTracks, inputParticleHypotheses,
      outputTrackParameters, outputSeeds, outputProtoTracks, trackingGeometry,
      magneticField, bFieldMin, initialSigmas, initialSigmaQoverPt,
      initialSigmaPtRel, initialVarInflation, noTimeVarInflation,
      particleHypothesis);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackParamsLookupEstimation, mex, "TrackParamsLookupEstimation",
      refLayers, bins, inputHits, inputParticles, trackLookupGridWriters);

  {
    using Alg = TrackFindingAlgorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, IAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackFindingAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config)
            .def_static("makeTrackFinderFunction",
                        [](std::shared_ptr<const Acts::TrackingGeometry>
                               trackingGeometry,
                           std::shared_ptr<const Acts::MagneticFieldProvider>
                               magneticField,
                           Acts::Logging::Level level) {
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
}

}  // namespace ActsPython
