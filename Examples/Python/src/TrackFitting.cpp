// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/TrackFitting/GsfFitterFunction.hpp"
#include "ActsExamples/TrackFitting/KalmanFitterFunction.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"

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

  {
    ACTS_PYTHON_DECLARE_ALGORITHM(
        ActsExamples::TrackFittingAlgorithm, mex, "TrackFittingAlgorithm",
        inputMeasurements, directNavigation, inputSourceLinks, inputProtoTracks,
        inputInitialTrackParameters, outputTrajectories, fit, trackingGeometry,
        pickTrack);

    mex.def(
        "makeKalmanFitterFunction",
        py::overload_cast<std::shared_ptr<const Acts::TrackingGeometry>,
                          std::shared_ptr<const Acts::MagneticFieldProvider>,
                          bool, bool, double, Acts::FreeToBoundCorrection>(
            &ActsExamples::makeKalmanFitterFunction),
        py::arg("trackingGeometry"), py::arg("magneticField"),
        py::arg("multipleScattering"), py::arg("energyLoss"),
        py::arg("reverseFilteringMomThreshold"),
        py::arg("freeToBoundCorrection"));

    py::enum_<Acts::FinalReductionMethod>(mex, "FinalReductionMethod")
        .value("mean", Acts::FinalReductionMethod::eMean)
        .value("maxWeight", Acts::FinalReductionMethod::eMaxWeight);

    mex.def(
        "makeGsfFitterFunction",
        py::overload_cast<std::shared_ptr<const Acts::TrackingGeometry>,
                          std::shared_ptr<const Acts::MagneticFieldProvider>, std::string, std::string,
                          std::size_t, Acts::FinalReductionMethod, bool, bool>(
            &ActsExamples::makeGsfFitterFunction),
        py::arg("trackingGeometry"), py::arg("magneticField"),
        py::arg("lowBetheHeitlerPath"), py::arg("highBetheHeitlerPath"),
        py::arg("maxComponents"), py::arg("finalReductionMethod"),
        py::arg("abortOnError"), py::arg("disableAllMaterialHandling"));
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
