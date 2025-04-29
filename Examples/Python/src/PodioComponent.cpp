// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/Podio/PodioReader.hpp"
#include "ActsExamples/Io/Podio/PodioWriter.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace Acts::Python;

PYBIND11_MODULE(ActsPythonBindingsPodio, m) {
  ACTS_PYTHON_DECLARE_READER(ActsExamples::PodioReader, m, "PodioReader",
                             inputPath, outputFrame, category);

  ACTS_PYTHON_DECLARE_WRITER(ActsExamples::PodioWriter, m, "PodioWriter",
                             inputFrame, outputPath, category, collections);
}
