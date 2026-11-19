// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/MeasurementFilterAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedParamsEstimator.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsLookupEstimation.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addSeeding(py::module& mex);

void addTrackFinding(py::module& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(MeasurementFilterAlgorithm, mex,
                                "MeasurementFilterAlgorithm", inputTracks,
                                inputMeasurementSubset, outputMeasurementSubset,
                                includeOutliers);

  addSeeding(mex);

  {
    auto e =
        py::class_<SeedParamsEstimator, std::shared_ptr<SeedParamsEstimator>>(
            mex, "SeedParamsEstimator");
    e.def_static("inverseRadiusPowerWeight",
                 &SeedParamsEstimator::inverseRadiusPowerWeight,
                 py::arg("exponent"));
    auto c =
        py::class_<SeedParamsEstimator::Config>(e, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(
        c, trackingGeometry, magneticField, bFieldMin, spacePointSelection,
        minTransverseDistance, geometricRefineIterations, spacePointWeight,
        initialSigmas, initialSigmaQoverPt, initialSigmaPtRel,
        initialVarInflation, noTimeVarInflation, particleHypothesis);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackParamsEstimationAlgorithm, mex, "TrackParamsEstimationAlgorithm",
      inputSeeds, inputProtoTracks, inputParticleHypotheses,
      outputTrackParameters, outputSeeds, outputProtoTracks,
      inputFreeParameters, seedParamsEstimation, trackingGeometry,
      magneticField);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackParamsLookupEstimation, mex, "TrackParamsLookupEstimation",
      refLayers, bins, inputHits, inputParticles, trackLookupGridWriters);

  {
    using Alg = TrackFindingAlgorithm;

    auto [alg, c] =
        declareAlgorithm<Alg, IAlgorithm>(mex, "TrackFindingAlgorithm");

    alg.def_static(
        "makeTrackFinderFunction",
        [](std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
           std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
           Acts::Logging::Level level) {
          return Alg::makeTrackFinderFunction(
              std::move(trackingGeometry), std::move(magneticField),
              *Acts::getDefaultLogger("TrackFinding", level));
        });

    alg.def_static(
        "makeBremTrackFinderFunction",
        [](std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
           std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
           Acts::Logging::Level level) {
          return Alg::makeBremTrackFinderFunction(
              std::move(trackingGeometry), std::move(magneticField),
              *Acts::getDefaultLogger("TrackFindingBrem", level));
        });

    py::class_<Alg::TrackFinderFunction,
               std::shared_ptr<Alg::TrackFinderFunction>>(
        alg, "TrackFinderFunction");

    py::class_<ActsExamples::TrackStateSelectionCuts>(mex,
                                                      "TrackStateSelectionCuts")
        .def(py::init<>())
        .def(py::init([](double chi2CutOffMeasurement, double chi2CutOffOutlier,
                         std::size_t maxBranchesPerSurface) {
               return ActsExamples::TrackStateSelectionCuts{
                   chi2CutOffMeasurement, chi2CutOffOutlier,
                   maxBranchesPerSurface};
             }),
             py::arg("chi2CutOffMeasurement") = 15.,
             py::arg("chi2CutOffOutlier") = 25.,
             py::arg("maxBranchesPerSurface") = 1)
        .def_readwrite(
            "chi2CutOffMeasurement",
            &ActsExamples::TrackStateSelectionCuts::chi2CutOffMeasurement)
        .def_readwrite(
            "chi2CutOffOutlier",
            &ActsExamples::TrackStateSelectionCuts::chi2CutOffOutlier)
        .def_readwrite(
            "maxBranchesPerSurface",
            &ActsExamples::TrackStateSelectionCuts::maxBranchesPerSurface);

    ACTS_PYTHON_STRUCT(
        c, inputMeasurements, inputInitialTrackParameters, inputSeeds,
        outputTracks, trackingGeometry, magneticField, findTracks,
        findTracksBrem, trackStateSelection, trackSelectorCfg, maxSteps, twoWay,
        twoWaySmoothing, seedParamsEstimation, reverseSearch, reverseSearchMinR,
        reverseSearchMinAbsZ, inwardOnlyMinR, inwardOnlyMinAbsZ,
        seedDeduplication, seedDeduplicatorCfg, stayOnSeed, pixelVolumeIds,
        stripVolumeIds, maxPixelHoles, maxStripHoles, trimTracks,
        recordMaterialStates, useJosephFormulation, constrainToVolumeIds,
        endOfWorldVolumeIds);
  }
}

}  // namespace ActsPython
