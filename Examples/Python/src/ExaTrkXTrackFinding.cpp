// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFindingOnnx.hpp"
#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFindingTorch.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFindingExaTrkX/ParameterFromTrajectoryAlgorithm.hpp"
#include "ActsExamples/TrackFindingExaTrkX/SourceLinkSelectorAlgorithm.hpp"
#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"
#include "ActsExamples/TrackFindingExaTrkX/TrackFindingFromPrototrackAlgorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addExaTrkXTrackFinding(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    using C = Acts::ExaTrkXTrackFindingBase;
    auto c = py::class_<C, std::shared_ptr<C>>(mex, "ExaTrkXTrackFindingBase");
  }

#ifdef ACTS_EXATRKX_TORCH_BACKEND
  {
    using Backend = Acts::ExaTrkXTrackFindingTorch;
    using Config = Backend::Config;

    auto backend =
        py::class_<Backend, Acts::ExaTrkXTrackFindingBase, std::shared_ptr<Backend>>(
            mex, "ExaTrkXTrackFindingTorch")
            .def(py::init<const Config&>(), py::arg("config"))
            .def_property_readonly("config", &Backend::config);

    auto c = py::class_<Config>(backend, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(modelDir);
    ACTS_PYTHON_MEMBER(spacepointFeatures);
    ACTS_PYTHON_MEMBER(embeddingDim);
    ACTS_PYTHON_MEMBER(rVal);
    ACTS_PYTHON_MEMBER(knnVal);
    ACTS_PYTHON_MEMBER(filterCut);
    ACTS_PYTHON_MEMBER(n_chunks);
    ACTS_PYTHON_MEMBER(edgeCut);
    ACTS_PYTHON_STRUCT_END();
  }
#endif

#ifdef ACTS_EXATRKX_ONNX_BACKEND
  {
    using Backend = Acts::ExaTrkXTrackFindingOnnx;
    using Config = Backend::Config;

    auto  =
        py::class_<Backend, Acts::ExaTrkXTrackFindingBase, std::shared_ptr<Backend>>(
            mex, "ExaTrkXTrackFindingOnnx")
            .def(py::init<const Config&>(), py::arg("config"))
            .def_property_readonly("config", &Backend::config);

    auto c = py::class_<Config>(, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(modelDir);
    ACTS_PYTHON_MEMBER(spacepointFeatures);
    ACTS_PYTHON_MEMBER(embeddingDim);
    ACTS_PYTHON_MEMBER(rVal);
    ACTS_PYTHON_MEMBER(knnVal);
    ACTS_PYTHON_MEMBER(filterCut);
    ACTS_PYTHON_STRUCT_END();
  }
#endif

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackFindingAlgorithmExaTrkX,
                                "TrackFindingAlgorithmExaTrkX",
                                inputSpacePoints, outputProtoTracks,
                                trackFinderML, rScale, phiScale, zScale)

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::SourceLinkSelectorAlgorithm,
                                "SourceLinkSelectorAlgorithm", inputSourceLinks,
                                outputSourceLinks, geometrySelection);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::ParameterFromTrajectoryAlgorithm,
                                "ParameterFromTrajectoryAlgorithm",
                                inputTrajectories, outputParamters)

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TrackFindingFromPrototrackAlgorithm,
      "TrackFindingFromPrototrackAlgorithm", inputTracks, inputMeasurements,
      inputSourceLinks, inputInitialTrackParameters, outputTrajectories,
      measurementSelectorCfg, trackingGeometry, magneticField)
}

}  // namespace Acts::Python
