// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/BoostTrackBuilding.hpp"
#include "Acts/Plugins/ExaTrkX/CugraphTrackBuilding.hpp"
#include "Acts/Plugins/ExaTrkX/OnnxEdgeClassifier.hpp"
#include "Acts/Plugins/ExaTrkX/OnnxMetricLearning.hpp"
#include "Acts/Plugins/ExaTrkX/TorchEdgeClassifier.hpp"
#include "Acts/Plugins/ExaTrkX/TorchMetricLearning.hpp"
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
    using C = Acts::GraphConstructionBase;
    auto c = py::class_<C, std::shared_ptr<C>>(mex, "GraphConstructionBase");
  }
  {
    using C = Acts::EdgeClassificationBase;
    auto c = py::class_<C, std::shared_ptr<C>>(mex, "EdgeClassificationBase");
  }
  {
    using C = Acts::TrackBuildingBase;
    auto c = py::class_<C, std::shared_ptr<C>>(mex, "TrackBuildingBase");
  }

#ifdef ACTS_EXATRKX_TORCH_BACKEND
  {
    using Alg = Acts::TorchMetricLearning;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, Acts::GraphConstructionBase, std::shared_ptr<Alg>>(
            mex, "TorchMetricLearning")
            .def(py::init<const Config&>(), py::arg("config"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(modelPath);
    ACTS_PYTHON_MEMBER(spacepointFeatures);
    ACTS_PYTHON_MEMBER(embeddingDim);
    ACTS_PYTHON_MEMBER(rVal);
    ACTS_PYTHON_MEMBER(knnVal);
    ACTS_PYTHON_STRUCT_END();
  }
  {
    using Alg = Acts::TorchEdgeClassifier;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, Acts::EdgeClassificationBase, std::shared_ptr<Alg>>(
            mex, "TorchEdgeClassifier")
            .def(py::init<const Config&>(), py::arg("config"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(modelPath);
    ACTS_PYTHON_MEMBER(cut);
    ACTS_PYTHON_MEMBER(nChunks);
    ACTS_PYTHON_STRUCT_END();
  }
  {
    using Alg = Acts::BoostTrackBuilding;

    auto alg = py::class_<Alg, Acts::TrackBuildingBase, std::shared_ptr<Alg>>(
                   mex, "BoostTrackBuilding")
                   .def(py::init<>());
  }
#endif

#ifdef ACTS_EXATRKX_ONNX_BACKEND
  {
    using Alg = Acts::OnnxMetricLearning;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, Acts::GraphConstructionBase, std::shared_ptr<Alg>>(
            mex, "OnnxMetricLearning")
            .def(py::init<const Config&>(), py::arg("config"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(modelPath);
    ACTS_PYTHON_MEMBER(spacepointFeatures);
    ACTS_PYTHON_MEMBER(embeddingDim);
    ACTS_PYTHON_MEMBER(rVal);
    ACTS_PYTHON_MEMBER(knnVal);
    ACTS_PYTHON_STRUCT_END();
  }
  {
    using Alg = Acts::OnnxEdgeClassifier;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, Acts::EdgeClassificationBase, std::shared_ptr<Alg>>(
            mex, "OnnxEdgeClassifier")
            .def(py::init<const Config&>(), py::arg("config"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(modelPath);
    ACTS_PYTHON_MEMBER(cut);
    ACTS_PYTHON_STRUCT_END();
  }
  {
    using Alg = Acts::CugraphTrackBuilding;

    auto alg = py::class_<Alg, Acts::TrackBuildingBase, std::shared_ptr<Alg>>(
                   mex, "CugraphTrackBuilding")
                   .def(py::init<>());
  }
#endif

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackFindingAlgorithmExaTrkX, mex,
                                "TrackFindingAlgorithmExaTrkX",
                                inputSpacePoints, outputProtoTracks,
                                graphConstructor, edgeClassifiers, trackBuilder,
                                rScale, phiScale, zScale);
}

}  // namespace Acts::Python
