// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFinding.hpp"
#include "ActsExamples/TrackFindingMLBased/TrackFindingMLBasedAlgorithm.hpp"

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
    using Alg = ActsExamples::TrackFindingMLBasedAlgorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackFindingMLBasedAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputSpacePoints);
    ACTS_PYTHON_MEMBER(outputProtoTracks);
    ACTS_PYTHON_MEMBER(exaTrkxConfig);
    ACTS_PYTHON_STRUCT_END();
  }

}

}  // namespace Acts::Python
