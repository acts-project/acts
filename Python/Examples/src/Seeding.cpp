// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsGraphBuilder.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"
#include "ActsExamples/TrackFinding/AdaptiveHoughTransformSeeder.hpp"
#include "ActsExamples/TrackFinding/GraphBasedSeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/GridTripletSeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/HoughTransformSeeder.hpp"
#include "ActsExamples/TrackFinding/MuonHoughSeeder.hpp"
#include "ActsExamples/TrackFinding/OrthogonalTripletSeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addSeeding(py::module& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(SpacePointMaker, mex, "SpacePointMaker",
                                inputMeasurements, outputSpacePoints,
                                trackingGeometry, geometrySelection,
                                stripGeometrySelection, stripVertex,
                                stripLengthTolerance, stripLengthGapTolerance);

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
      maxQualitySeedsPerSpMConf, useDeltaRinsteadOfTopRadius, useTimeDoubletCut,
      doubletTimeCutNSigma, useTimeTripletCut, tripletTimeChi2Max, useExtraCuts,
      inputVertices, vertexZNSigma, vertexZMargin);

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

  {
    using Config = Acts::Experimental::GraphBasedTrackSeeder::Config;
    auto c =
        py::class_<Config>(mex, "GraphBasedSeedingConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, nMaxPhiSlice, maxSeedSplitEta,
                       postFitCurvatureFactor);
    patchKwargsConstructor(c);
  }

  {
    using Config = Acts::Experimental::GbtsGraphBuilder::Config;
    auto c =
        py::class_<Config>(mex, "GbtsGraphBuilderConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, minPt, maxEdgesPerSP, minDeltaRadius, doubletFilterRZ,
                       matchBeforeCreate, validateTriplets, calibrateStrips,
                       maxOuterRadius, z0HistogramMaxBarrelOrder,
                       useSumRatioCut, sumRatioCutMax);
    patchKwargsConstructor(c);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      GraphBasedSeedingAlgorithm, mex, "GraphBasedSeedingAlgorithm",
      inputSpacePoints, outputSeeds, seedFinderConfig, graphConfig,
      useStripConnections, useAdaptiveCuts, layerMappingFile,
      connectorInputFile, lutInputFile, collisionRegionMin, collisionRegionMax,
      trackingGeometry, fillModuleCsv, inputClusters, bFieldInZ, useLayerTilt,
      magneticField, outputFreeParameters);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      HoughTransformSeeder, mex, "HoughTransformSeeder", inputSpacePoints,
      outputProtoTracks, trackingGeometry, geometrySelection, inputMeasurements,
      subRegions, nLayers, xMin, xMax, yMin, yMax, houghHistSize_x,
      houghHistSize_y, hitExtend_x, threshold, localMaxWindowSize, kA);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      AdaptiveHoughTransformSeeder, mex, "AdaptiveHoughTransformSeeder",
      inputSpacePoints, outputSeeds, trackingGeometry, qOverPtMin,
      qOverPtMinBinSize, phiMinBinSize, threshold, noiseThreshold, deduplicate,
      inverseA, doSecondPhase, zRange, cotThetaRange, cotThetaMinBinSize,
      zMinBinSize);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      MuonHoughSeeder, mex, "MuonHoughSeeder", inTruthSegments, inSpacePoints,
      outHoughMax, nBinsTanTheta, nBinsY0, nBinsTanPhi, nBinsX0,
      dumpVisualization, extendWithPhi, visualizationFunction);
}

}  // namespace ActsPython
