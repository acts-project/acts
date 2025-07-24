// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Printers/ParticlesPrinter.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include "ActsPython/Utilities/Patchers.hpp"

#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsPython;
using namespace ActsExamples;
using namespace Acts;

namespace ActsPython {

/// This adds the printer algorithms to the examples module
/// @param mex the examples module
void addPrinterAlgorithms(py::module_& mex) {
  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::ParticlesPrinter, mex,
                                "ParticlesPrinter", inputParticles);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackParametersPrinter, mex,
                                "TrackParametersPrinter", inputTrackParameters);
}
}  // namespace ActsPython
