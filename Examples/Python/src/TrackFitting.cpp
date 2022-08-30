// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
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

  {
    using Alg = ActsExamples::SurfaceSortingAlgorithm;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "SurfaceSortingAlgorithm")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputProtoTracks);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(inputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(outputProtoTracks);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::TrackFittingAlgorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackFittingAlgorithm")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config)
            .def_static("makeKalmanFitterFunction",
                        py::overload_cast<
                            std::shared_ptr<const Acts::TrackingGeometry>,
                            std::shared_ptr<const Acts::MagneticFieldProvider>,
                            bool, bool, double, Acts::FreeToBoundCorrection>(
                            &Alg::makeKalmanFitterFunction),
                        py::arg("trackingGeometry"), py::arg("magneticField"),
                        py::arg("multipleScattering"), py::arg("energyLoss"),
                        py::arg("reverseFilteringMomThreshold"),
                        py::arg("freeToBoundCorrection"))
            .def_static("makeKalmanFitterFunction",
                        py::overload_cast<
                            std::shared_ptr<const Acts::MagneticFieldProvider>,
                            bool, bool, double, Acts::FreeToBoundCorrection>(
                            &Alg::makeKalmanFitterFunction),
                        py::arg("magneticField"), py::arg("multipleScattering"),
                        py::arg("energyLoss"),
                        py::arg("reverseFilteringMomThreshold"),
                        py::arg("freeToBoundCorrection"))
            .def_static(
                "makeGsfFitterFunction",
                py::overload_cast<
                    std::shared_ptr<const Acts::TrackingGeometry>,
                    std::shared_ptr<const Acts::MagneticFieldProvider>,
                    std::size_t, bool, bool>(&Alg::makeGsfFitterFunction),
                py::arg("trackingGeometry"), py::arg("magneticField"),
                py::arg("maxComponents"), py::arg("abortOnError"),
                py::arg("disableAllMaterialHandling"))
            .def_static(
                "makeGsfFitterFunction",
                py::overload_cast<
                    std::shared_ptr<const Acts::MagneticFieldProvider>,
                    std::size_t, bool, bool>(&Alg::makeGsfFitterFunction),
                py::arg("magneticField"), py::arg("maxComponents"),
                py::arg("abortOnError"), py::arg("disableAllMaterialHandling"));

    py::class_<TrackFittingAlgorithm::TrackFitterFunction,
               std::shared_ptr<TrackFittingAlgorithm::TrackFitterFunction>>(
        alg, "TrackFitterFunction");

    py::class_<
        TrackFittingAlgorithm::DirectedTrackFitterFunction,
        std::shared_ptr<TrackFittingAlgorithm::DirectedTrackFitterFunction>>(
        alg, "DirectedTrackFitterFunction");

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(directNavigation);
    ACTS_PYTHON_MEMBER(inputSourceLinks);
    ACTS_PYTHON_MEMBER(inputProtoTracks);
    ACTS_PYTHON_MEMBER(inputInitialTrackParameters);
    ACTS_PYTHON_MEMBER(outputTrajectories);
    ACTS_PYTHON_MEMBER(fit);
    ACTS_PYTHON_MEMBER(dFit);
    ACTS_PYTHON_MEMBER(trackingGeometry);
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
