// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Seeding/SeedFinderOrthogonalConfig.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "ActsExamples/TrackFinding/AmbiguityResolutionAlgorithm.hpp"
#include "ActsExamples/TrackFinding/HoughTransformSeeder.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedingOrthogonalAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrajectoriesToPrototracks.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addTrackFinding(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::SpacePointMaker, mex,
                                "SpacePointMaker", inputSourceLinks,
                                inputMeasurements, outputSpacePoints,
                                trackingGeometry, geometrySelection);

  {
    using Config = Acts::SeedFilterConfig;
    auto c = py::class_<Config>(m, "SeedFilterConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(deltaInvHelixDiameter);
    ACTS_PYTHON_MEMBER(impactWeightFactor);
    ACTS_PYTHON_MEMBER(zOriginWeightFactor);
    ACTS_PYTHON_MEMBER(compatSeedWeight);
    ACTS_PYTHON_MEMBER(deltaRMin);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpM);
    ACTS_PYTHON_MEMBER(compatSeedLimit);
    ACTS_PYTHON_MEMBER(curvatureSortingInFilter);
    ACTS_PYTHON_MEMBER(seedConfirmation);
    ACTS_PYTHON_MEMBER(centralSeedConfirmationRange);
    ACTS_PYTHON_MEMBER(forwardSeedConfirmationRange);
    ACTS_PYTHON_MEMBER(useDeltaRorTopRadius);
    ACTS_PYTHON_MEMBER(seedWeightIncrement);
    ACTS_PYTHON_MEMBER(numSeedIncrement);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpMConf);
    ACTS_PYTHON_MEMBER(maxQualitySeedsPerSpMConf);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }

  {
    using Config = Acts::SeedFinderConfig<SimSpacePoint>;
    auto c = py::class_<Config>(m, "SeedFinderConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(cotThetaMax);
    ACTS_PYTHON_MEMBER(deltaRMin);
    ACTS_PYTHON_MEMBER(deltaRMax);
    ACTS_PYTHON_MEMBER(deltaRMinBottomSP);
    ACTS_PYTHON_MEMBER(deltaRMaxBottomSP);
    ACTS_PYTHON_MEMBER(deltaRMinTopSP);
    ACTS_PYTHON_MEMBER(deltaRMaxTopSP);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_MEMBER(sigmaScattering);
    ACTS_PYTHON_MEMBER(maxPtScattering);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpM);
    ACTS_PYTHON_MEMBER(collisionRegionMin);
    ACTS_PYTHON_MEMBER(collisionRegionMax);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(rMax);
    ACTS_PYTHON_MEMBER(rMin);
    ACTS_PYTHON_MEMBER(radLengthPerSeed);
    ACTS_PYTHON_MEMBER(zAlign);
    ACTS_PYTHON_MEMBER(rAlign);
    ACTS_PYTHON_MEMBER(sigmaError);
    ACTS_PYTHON_MEMBER(maxBlockSize);
    ACTS_PYTHON_MEMBER(nTrplPerSpBLimit);
    ACTS_PYTHON_MEMBER(nAvgTrplPerSpBLimit);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_MEMBER(deltaZMax);
    ACTS_PYTHON_MEMBER(zBinEdges);
    ACTS_PYTHON_MEMBER(skipPreviousTopSP);
    ACTS_PYTHON_MEMBER(interactionPointCut);
    ACTS_PYTHON_MEMBER(zBinsCustomLooping);
    ACTS_PYTHON_MEMBER(useVariableMiddleSPRange);
    ACTS_PYTHON_MEMBER(deltaRMiddleMinSPRange);
    ACTS_PYTHON_MEMBER(deltaRMiddleMaxSPRange);
    ACTS_PYTHON_MEMBER(rRangeMiddleSP);
    ACTS_PYTHON_MEMBER(rMinMiddle);
    ACTS_PYTHON_MEMBER(rMaxMiddle);
    ACTS_PYTHON_MEMBER(binSizeR);
    ACTS_PYTHON_MEMBER(forceRadialSorting);
    ACTS_PYTHON_MEMBER(seedConfirmation);
    ACTS_PYTHON_MEMBER(centralSeedConfirmationRange);
    ACTS_PYTHON_MEMBER(forwardSeedConfirmationRange);
    ACTS_PYTHON_MEMBER(arithmeticAverageCotTheta);
    ACTS_PYTHON_MEMBER(useDetailedDoubleMeasurementInfo);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }
  {
    using seedOptions = Acts::SeedFinderOptions;
    auto c = py::class_<seedOptions>(m, "SeedFinderOptions").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, seedOptions);
    ACTS_PYTHON_MEMBER(beamPos);
    ACTS_PYTHON_MEMBER(bFieldInZ);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }
  {
    using Config = Acts::SeedFinderOrthogonalConfig<SimSpacePoint>;
    auto c =
        py::class_<Config>(m, "SeedFinderOrthogonalConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(cotThetaMax);
    ACTS_PYTHON_MEMBER(deltaRMinBottomSP);
    ACTS_PYTHON_MEMBER(deltaRMaxBottomSP);
    ACTS_PYTHON_MEMBER(deltaRMinTopSP);
    ACTS_PYTHON_MEMBER(deltaRMaxTopSP);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_MEMBER(deltaZMax);
    ACTS_PYTHON_MEMBER(sigmaScattering);
    ACTS_PYTHON_MEMBER(maxPtScattering);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpM);
    ACTS_PYTHON_MEMBER(collisionRegionMin);
    ACTS_PYTHON_MEMBER(collisionRegionMax);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(rMax);
    ACTS_PYTHON_MEMBER(rMin);
    ACTS_PYTHON_MEMBER(radLengthPerSeed);
    ACTS_PYTHON_MEMBER(deltaZMax);
    ACTS_PYTHON_MEMBER(skipPreviousTopSP);
    ACTS_PYTHON_MEMBER(interactionPointCut);
    ACTS_PYTHON_MEMBER(deltaPhiMax);
    ACTS_PYTHON_MEMBER(highland);
    ACTS_PYTHON_MEMBER(maxScatteringAngle2);
    ACTS_PYTHON_MEMBER(useVariableMiddleSPRange);
    ACTS_PYTHON_MEMBER(deltaRMiddleMinSPRange);
    ACTS_PYTHON_MEMBER(deltaRMiddleMaxSPRange);
    ACTS_PYTHON_MEMBER(rRangeMiddleSP);
    ACTS_PYTHON_MEMBER(rMinMiddle);
    ACTS_PYTHON_MEMBER(rMaxMiddle);
    ACTS_PYTHON_MEMBER(seedConfirmation);
    ACTS_PYTHON_MEMBER(centralSeedConfirmationRange);
    ACTS_PYTHON_MEMBER(forwardSeedConfirmationRange);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }

  {
    using seedConf = Acts::SeedConfirmationRangeConfig;
    auto c = py::class_<seedConf>(m, "SeedConfirmationRangeConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, seedConf);
    ACTS_PYTHON_MEMBER(zMinSeedConf);
    ACTS_PYTHON_MEMBER(zMaxSeedConf);
    ACTS_PYTHON_MEMBER(rMaxSeedConf);
    ACTS_PYTHON_MEMBER(nTopForLargeR);
    ACTS_PYTHON_MEMBER(nTopForSmallR);
    ACTS_PYTHON_MEMBER(seedConfMinBottomRadius);
    ACTS_PYTHON_MEMBER(seedConfMaxZOrigin);
    ACTS_PYTHON_MEMBER(minImpactSeedConf);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }

  {
    using Config = Acts::SpacePointGridConfig;
    auto c = py::class_<Config>(m, "SpacePointGridConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(rMax);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(deltaRMax);
    ACTS_PYTHON_MEMBER(cotThetaMax);
    ACTS_PYTHON_MEMBER(phiBinDeflectionCoverage);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_MEMBER(zBinEdges);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }
  {
    using Options = Acts::SpacePointGridOptions;
    auto c = py::class_<Options>(m, "SpacePointGridOptions").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Options);
    ACTS_PYTHON_MEMBER(bFieldInZ);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::SeedingAlgorithm, mex, "SeedingAlgorithm", inputSpacePoints,
      outputSeeds, outputProtoTracks, seedFilterConfig, seedFinderConfig,
      seedFinderOptions, gridConfig, gridOptions, allowSeparateRMax,
      zBinNeighborsTop, zBinNeighborsBottom, numPhiNeighbors);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::SeedingOrthogonalAlgorithm, mex,
      "SeedingOrthogonalAlgorithm", inputSpacePoints, outputSeeds,
      outputProtoTracks, seedFilterConfig, seedFinderConfig, seedFinderOptions);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::HoughTransformSeeder, mex, "HoughTransformSeeder",
      inputSpacePoints, outputSeeds, outputProtoTracks, inputSourceLinks,
      trackingGeometry, geometrySelection, inputMeasurements, subRegions,
      nLayers, xMin, xMax, yMin, yMax, houghHistSize_x, houghHistSize_y,
      hitExtend_x, threshold, localMaxWindowSize, kA);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TrackParamsEstimationAlgorithm, mex,
      "TrackParamsEstimationAlgorithm", inputSeeds, inputSpacePoints,
      inputProtoTracks, inputSourceLinks, outputTrackParameters,
      outputProtoTracks, trackingGeometry, magneticField, deltaRMin, deltaRMax,
      bFieldMin, sigmaLoc0, sigmaLoc1, sigmaPhi, sigmaTheta, sigmaQOverP,
      sigmaT0, initialVarInflation);

  {
    using Alg = ActsExamples::TrackFindingAlgorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackFindingAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config)
            .def_static("makeTrackFinderFunction",
                        &Alg::makeTrackFinderFunction);

    py::class_<Alg::TrackFinderFunction,
               std::shared_ptr<Alg::TrackFinderFunction>>(
        alg, "TrackFinderFunction");

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(inputSourceLinks);
    ACTS_PYTHON_MEMBER(inputInitialTrackParameters);
    ACTS_PYTHON_MEMBER(outputTrajectories);
    ACTS_PYTHON_MEMBER(findTracks);
    ACTS_PYTHON_MEMBER(measurementSelectorCfg);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::TrajectoriesToPrototracks;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrajectoriesToPrototracks")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputTrajectories);
    ACTS_PYTHON_MEMBER(outputPrototracks);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    auto constructor = [](const std::vector<
                           std::pair<GeometryIdentifier,
                                     std::tuple<std::vector<double>,
                                                std::vector<double>,
                                                std::vector<size_t>>>>& input) {
      std::vector<std::pair<GeometryIdentifier, MeasurementSelectorCuts>>
          converted;
      converted.reserve(input.size());
      for (const auto& [id, cuts] : input) {
        const auto& [bins, chi2, num] = cuts;
        converted.emplace_back(id, MeasurementSelectorCuts{bins, chi2, num});
      }
      return std::make_unique<MeasurementSelector::Config>(converted);
    };

    py::class_<MeasurementSelectorCuts>(m, "MeasurementSelectorCuts")
        .def(py::init<>())
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<size_t>>())
        .def_readwrite("etaBins", &MeasurementSelectorCuts::etaBins)
        .def_readwrite("chi2CutOff", &MeasurementSelectorCuts::chi2CutOff)
        .def_readwrite("numMeasurementsCutOff",
                       &MeasurementSelectorCuts::numMeasurementsCutOff);

    auto ms = py::class_<MeasurementSelector>(m, "MeasurementSelector");
    auto c =
        py::class_<MeasurementSelector::Config>(ms, "Config")
            .def(py::init<std::vector<
                     std::pair<GeometryIdentifier, MeasurementSelectorCuts>>>())
            .def(py::init(constructor));
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::AmbiguityResolutionAlgorithm, mex,
                                "AmbiguityResolutionAlgorithm",
                                inputSourceLinks, inputTrajectories,
                                outputTrajectories, maximumSharedHits);
}

}  // namespace Acts::Python
