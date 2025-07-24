// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Hashing/HashingAlgorithmConfig.hpp"
#include "Acts/Plugins/Hashing/HashingTrainingConfig.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include "ActsPython/Utilities/Patchers.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace Acts;

namespace ActsPython {

/// This adds the hashing algorithms to the plugins
/// @param p the plugins module
void addHashing(py::module_& p) {

  auto hashing = p.def_submodule("hashing");

  {
    using Config = Acts::HashingAlgorithmConfig;
    auto c = py::class_<Config>(hashing, "HashingAlgorithmConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, bucketSize, zBins, phiBins);
    patchKwargsConstructor(c);
  }

  {
    using Config = Acts::HashingTrainingConfig;
    auto c = py::class_<Config>(hashing, "HashingTrainingConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, annoySeed, f);
    patchKwargsConstructor(c);
  }
}

}  // namespace ActsPython