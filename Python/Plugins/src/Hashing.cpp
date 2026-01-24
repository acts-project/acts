// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Hashing/HashingAlgorithm.hpp"
#include "ActsPlugins/Hashing/HashingTraining.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsHashing, hashing) {
  using namespace ActsPlugins;
  using namespace ActsPython;

  {
    using Config = HashingAlgorithm::Config;
    auto c =
        py::class_<Config>(hashing, "HashingAlgorithmConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, bucketSize, zBins, phiBins);
    patchKwargsConstructor(c);
  }

  {
    using Config = HashingTraining::Config;
    auto c =
        py::class_<Config>(hashing, "HashingTrainingConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, annoySeed, f);
    patchKwargsConstructor(c);
  }
}
