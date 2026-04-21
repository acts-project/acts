// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/AlignmentMillePede/ActsSolverFromMille.hpp"
#include "ActsExamples/AlignmentMillePede/MillePedeAlignmentSandbox.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsAlignmentMillePede, m) {
  ACTS_PYTHON_DECLARE_ALGORITHM(MillePedeAlignmentSandbox, m,
                                "MillePedeAlignmentSandbox", milleOutput,
                                inputMeasurements, inputTracks,
                                trackingGeometry, magneticField, fixModules);
  ACTS_PYTHON_DECLARE_ALGORITHM(ActsSolverFromMille, m, "ActsSolverFromMille",
                                milleInput, txtOutput, trackingGeometry,
                                magneticField, fixModules);
}
