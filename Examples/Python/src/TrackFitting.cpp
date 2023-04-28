// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/TrackFitting/RefittingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TrackFittingChi2/TrackFittingChi2Algorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TrackFittingAlgorithm, mex, "TrackFittingAlgorithm",
      inputMeasurements, inputSourceLinks, inputProtoTracks,
      inputInitialTrackParameters, outputTracks, fit, pickTrack, calibrator);

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

    py::enum_<Acts::FinalReductionMethod>(mex, "FinalReductionMethod")
        .value("mean", Acts::FinalReductionMethod::eMean)
        .value("maxWeight", Acts::FinalReductionMethod::eMaxWeight);

    py::class_<ActsExamples::BetheHeitlerApprox>(mex, "AtlasBetheHeitlerApprox")
        .def_static("loadFromFiles",
                    &ActsExamples::BetheHeitlerApprox::loadFromFiles,
                    py::arg("lowParametersPath"), py::arg("lowParametersPath"))
        .def_static("makeDefault", []() {
          return Acts::Experimental::makeDefaultBetheHeitlerApprox();
        });

    mex.def(
        "makeGsfFitterFunction",
        [](std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
           std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
           BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
           double weightCutoff, Acts::FinalReductionMethod finalReductionMethod,
           bool abortOnError, bool disableAllMaterialHandling,
           Logging::Level level) {
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
  }

  {
    using Alg = ActsExamples::TrackFittingChi2Algorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, IAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackFittingChi2Algorithm")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config)
            .def_static("makeTrackFitterChi2Function",
                        py::overload_cast<
                            std::shared_ptr<const Acts::TrackingGeometry>,
                            std::shared_ptr<const Acts::MagneticFieldProvider>>(
                            &Alg::makeTrackFitterChi2Function));

    py::class_<
        TrackFittingChi2Algorithm::TrackFitterChi2Function,
        std::shared_ptr<TrackFittingChi2Algorithm::TrackFitterChi2Function>>(
        alg, "TrackFitterChi2Function");

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(inputSourceLinks);
    ACTS_PYTHON_MEMBER(inputProtoTracks);
    ACTS_PYTHON_MEMBER(inputInitialTrackParameters);
    ACTS_PYTHON_MEMBER(outputTracks);
    ACTS_PYTHON_MEMBER(nUpdates);
    ACTS_PYTHON_MEMBER(fit);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(multipleScattering);
    ACTS_PYTHON_MEMBER(energyLoss);
    ACTS_PYTHON_MEMBER(pickTrack);
    ACTS_PYTHON_STRUCT_END();
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
