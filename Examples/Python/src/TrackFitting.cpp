// This file is part of the Acts project.
//
// Copyright (C) 2021-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/Utilities/GaussianMixtureReduction.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/ScalingCalibrator.hpp"
#include "ActsExamples/TrackFitting/RefittingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"

#include <cstddef>
#include <memory>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace Acts {
class MagneticFieldProvider;
class TrackingGeometry;
}  // namespace Acts
namespace ActsExamples {
class IAlgorithm;
}  // namespace ActsExamples

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addTrackFitting(Context& ctx) {
  auto mex = ctx.get("examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::SurfaceSortingAlgorithm, mex,
                                "SurfaceSortingAlgorithm", inputProtoTracks,
                                inputSimHits, inputMeasurementSimHitsMap,
                                outputProtoTracks);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackFittingAlgorithm, mex,
                                "TrackFittingAlgorithm", inputMeasurements,
                                inputSourceLinks, inputProtoTracks,
                                inputInitialTrackParameters, inputClusters,
                                outputTracks, fit, pickTrack, calibrator);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::RefittingAlgorithm, mex,
                                "RefittingAlgorithm", inputTracks, outputTracks,
                                fit, pickTrack);

  {
    py::class_<TrackFitterFunction, std::shared_ptr<TrackFitterFunction>>(
        mex, "TrackFitterFunction");

    mex.def(
        "makeKalmanFitterFunction",
        [](std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
           std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
           bool multipleScattering, bool energyLoss,
           double reverseFilteringMomThreshold,
           Acts::FreeToBoundCorrection freeToBoundCorrection,
           Logging::Level level) {
          return ActsExamples::makeKalmanFitterFunction(
              trackingGeometry, magneticField, multipleScattering, energyLoss,
              reverseFilteringMomThreshold, freeToBoundCorrection,
              *Acts::getDefaultLogger("Kalman", level));
        },
        py::arg("trackingGeometry"), py::arg("magneticField"),
        py::arg("multipleScattering"), py::arg("energyLoss"),
        py::arg("reverseFilteringMomThreshold"),
        py::arg("freeToBoundCorrection"), py::arg("level"));

    py::class_<MeasurementCalibrator, std::shared_ptr<MeasurementCalibrator>>(
        mex, "MeasurementCalibrator");

    mex.def("makePassThroughCalibrator",
            []() -> std::shared_ptr<MeasurementCalibrator> {
              return std::make_shared<PassThroughCalibrator>();
            });

    mex.def(
        "makeScalingCalibrator",
        [](const char* path) -> std::shared_ptr<MeasurementCalibrator> {
          return std::make_shared<ActsExamples::ScalingCalibrator>(path);
        },
        py::arg("path"));

    py::enum_<Acts::MixtureReductionMethod>(mex, "FinalReductionMethod")
        .value("mean", Acts::MixtureReductionMethod::eMean)
        .value("maxWeight", Acts::MixtureReductionMethod::eMaxWeight);

    py::class_<ActsExamples::BetheHeitlerApprox>(mex, "AtlasBetheHeitlerApprox")
        .def_static("loadFromFiles",
                    &ActsExamples::BetheHeitlerApprox::loadFromFiles,
                    py::arg("lowParametersPath"), py::arg("lowParametersPath"))
        .def_static("makeDefault",
                    []() { return Acts::makeDefaultBetheHeitlerApprox(); });

    mex.def(
        "makeGsfFitterFunction",
        [](std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
           std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
           BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
           double weightCutoff,
           Acts::MixtureReductionMethod finalReductionMethod, bool abortOnError,
           bool disableAllMaterialHandling, Logging::Level level) {
          return ActsExamples::makeGsfFitterFunction(
              trackingGeometry, magneticField, betheHeitlerApprox,
              maxComponents, weightCutoff, finalReductionMethod, abortOnError,
              disableAllMaterialHandling,
              *Acts::getDefaultLogger("GSFFunc", level));
        },
        py::arg("trackingGeometry"), py::arg("magneticField"),
        py::arg("betheHeitlerApprox"), py::arg("maxComponents"),
        py::arg("weightCutoff"), py::arg("finalReductionMethod"),
        py::arg("abortOnError"), py::arg("disableAllMaterialHandling"),
        py::arg("level"));

    mex.def(
        "makeGlobalChiSquareFitterFunction",
        [](std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
           std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
           bool multipleScattering, bool energyLoss,
           Acts::FreeToBoundCorrection freeToBoundCorrection,
           Logging::Level level) {
          return ActsExamples::makeGlobalChiSquareFitterFunction(
              trackingGeometry, magneticField, multipleScattering, energyLoss,
              freeToBoundCorrection, *Acts::getDefaultLogger("Gx2f", level));
        },
        py::arg("trackingGeometry"), py::arg("magneticField"),
        py::arg("multipleScattering"), py::arg("energyLoss"),
        py::arg("freeToBoundCorrection"), py::arg("level"));

    // TODO add other important parameters like nUpdates
    // TODO add also in trackfitterfunction
  }

  {
    py::class_<FreeToBoundCorrection>(mex, "FreeToBoundCorrection")
        .def(py::init<>())
        .def(py::init<bool>(), py::arg("apply") = false)
        .def(py::init<bool, double, double>(), py::arg("apply") = false,
             py::arg("alpha") = 0.1, py::arg("beta") = 2);
  }
}

}  // namespace Acts::Python
