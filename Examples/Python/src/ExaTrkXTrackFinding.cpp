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
#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"

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
    using Alg = Acts::ExaTrkXTrackFindingTorch;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, Acts::ExaTrkXTrackFindingBase, std::shared_ptr<Alg>>(
            mex, "ExaTrkXTrackFindingTorch")
            .def(py::init<const Config&>(), py::arg("config"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
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
    using Alg = Acts::ExaTrkXTrackFindingOnnx;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, Acts::ExaTrkXTrackFindingBase, std::shared_ptr<Alg>>(
            mex, "ExaTrkXTrackFindingOnnx")
            .def(py::init<const Config&>(), py::arg("config"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
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

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackFindingAlgorithmExaTrkX, mex,
                                "TrackFindingAlgorithmExaTrkX",
                                inputSpacePoints, outputProtoTracks,
                                trackFinderML, rScale, phiScale, zScale);
}

}  // namespace Acts::Python
