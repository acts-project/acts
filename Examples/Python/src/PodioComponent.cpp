// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Podio/PodioReader.hpp"
#include "ActsExamples/Io/Podio/PodioWriter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsPythonBindingsPodio, m) {
  ACTS_PYTHON_DECLARE_READER(PodioReader, m, "PodioReader", inputPath,
                             outputFrame, category);

  ACTS_PYTHON_DECLARE_WRITER(PodioWriter, m, "PodioWriter", inputFrame,
                             outputPath, category, collections);
}
