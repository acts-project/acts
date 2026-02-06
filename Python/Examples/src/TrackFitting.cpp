// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/ScalingCalibrator.hpp"
#include "ActsExamples/TrackFitting/RefittingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <cstddef>
#include <memory>
#include <utility>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;
using namespace py::literals;

namespace ActsPython {

void addTrackFitting(py::module& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(
      TrackFittingAlgorithm, mex, "TrackFittingAlgorithm", inputMeasurements,
      inputProtoTracks, inputInitialTrackParameters, inputClusters,
      outputTracks, fit, pickTrack, calibrator);

  ACTS_PYTHON_DECLARE_ALGORITHM(RefittingAlgorithm, mex, "RefittingAlgorithm",
                                inputTracks, outputTracks, fit, pickTrack,
                                initialVarInflation);

  {
    py::class_<TrackFitterFunction, std::shared_ptr<TrackFitterFunction>>(
        mex, "TrackFitterFunction");

    mex.def(
        "makeKalmanFitterFunction",
        [](std::shared_ptr<const TrackingGeometry> trackingGeometry,
           std::shared_ptr<const MagneticFieldProvider> magneticField,
           bool multipleScattering, bool energyLoss,
           double reverseFilteringMomThreshold,
           double reverseFilteringCovarianceScaling,
           FreeToBoundCorrection freeToBoundCorrection, double chi2Cut,
           Logging::Level level) {
          return makeKalmanFitterFunction(
              std::move(trackingGeometry), std::move(magneticField),
              multipleScattering, energyLoss, reverseFilteringMomThreshold,
              reverseFilteringCovarianceScaling, freeToBoundCorrection, chi2Cut,
              *getDefaultLogger("Kalman", level));
        },
        "trackingGeometry"_a, "magneticField"_a, "multipleScattering"_a,
        "energyLoss"_a, "reverseFilteringMomThreshold"_a,
        "reverseFilteringCovarianceScaling"_a, "freeToBoundCorrection"_a,
        "chi2Cut"_a, "level"_a);

    py::class_<MeasurementCalibrator, std::shared_ptr<MeasurementCalibrator>>(
        mex, "MeasurementCalibrator");

    mex.def("makePassThroughCalibrator",
            []() -> std::shared_ptr<MeasurementCalibrator> {
              return std::make_shared<PassThroughCalibrator>();
            });

    mex.def(
        "makeScalingCalibrator",
        [](const char* path) -> std::shared_ptr<MeasurementCalibrator> {
          return std::make_shared<ScalingCalibrator>(path);
        },
        py::arg("path"));

    py::enum_<ComponentMergeMethod>(mex, "ComponentMergeMethod")
        .value("mean", ComponentMergeMethod::eMean)
        .value("maxWeight", ComponentMergeMethod::eMaxWeight);

    py::enum_<MixtureReductionAlgorithm>(mex, "MixtureReductionAlgorithm")
        .value("weightCut", MixtureReductionAlgorithm::weightCut)
        .value("KLDistance", MixtureReductionAlgorithm::KLDistance);

    py::class_<BetheHeitlerApprox, std::shared_ptr<BetheHeitlerApprox>>(
        mex, "BetheHeitlerApprox");
    py::class_<AtlasBetheHeitlerApprox, BetheHeitlerApprox,
               std::shared_ptr<AtlasBetheHeitlerApprox>>(
        mex, "AtlasBetheHeitlerApprox")
        .def_static("loadFromFiles", &AtlasBetheHeitlerApprox::loadFromFiles,
                    "lowParametersPath"_a, "highParametersPath"_a, "lowLimit"_a,
                    "highLimit"_a, "clampToRange"_a, "noChangeLimit"_a,
                    "singleGaussianLimit"_a)
        .def_static(
            "makeDefault",
            [](bool clampToRange) {
              return makeDefaultBetheHeitlerApprox(clampToRange);
            },
            "clampToRange"_a);

    mex.def(
        "makeGsfFitterFunction",
        [](std::shared_ptr<const TrackingGeometry> trackingGeometry,
           std::shared_ptr<const MagneticFieldProvider> magneticField,
           const std::shared_ptr<const BetheHeitlerApprox>& betheHeitlerApprox,
           std::size_t maxComponents, double weightCutoff,
           ComponentMergeMethod componentMergeMethod,
           MixtureReductionAlgorithm mixtureReductionAlgorithm,
           double reverseFilteringCovarianceScaling, Logging::Level level) {
          return makeGsfFitterFunction(
              std::move(trackingGeometry), std::move(magneticField),
              betheHeitlerApprox, maxComponents, weightCutoff,
              componentMergeMethod, mixtureReductionAlgorithm,
              reverseFilteringCovarianceScaling,
              *getDefaultLogger("GSFFunc", level));
        },
        "trackingGeometry"_a, "magneticField"_a, "betheHeitlerApprox"_a,
        "maxComponents"_a, "weightCutoff"_a, "componentMergeMethod"_a,
        "mixtureReductionAlgorithm"_a, "reverseFilteringCovarianceScaling"_a,
        "level"_a);

    mex.def(
        "makeGlobalChiSquareFitterFunction",
        [](std::shared_ptr<const TrackingGeometry> trackingGeometry,
           std::shared_ptr<const MagneticFieldProvider> magneticField,
           bool multipleScattering, bool energyLoss,
           FreeToBoundCorrection freeToBoundCorrection, std::size_t nUpdateMax,
           double relChi2changeCutOff, Logging::Level level) {
          return makeGlobalChiSquareFitterFunction(
              std::move(trackingGeometry), std::move(magneticField),
              multipleScattering, energyLoss, freeToBoundCorrection, nUpdateMax,
              relChi2changeCutOff, *getDefaultLogger("Gx2f", level));
        },
        py::arg("trackingGeometry"), py::arg("magneticField"),
        py::arg("multipleScattering"), py::arg("energyLoss"),
        py::arg("freeToBoundCorrection"), py::arg("nUpdateMax"),
        py::arg("relChi2changeCutOff"), py::arg("level"));
  }

  {
    py::class_<FreeToBoundCorrection>(mex, "FreeToBoundCorrection")
        .def(py::init<>())
        .def(py::init<bool>(), py::arg("apply") = false)
        .def(py::init<bool, double, double>(), py::arg("apply") = false,
             py::arg("alpha") = 0.1, py::arg("beta") = 2);
  }
}

}  // namespace ActsPython
