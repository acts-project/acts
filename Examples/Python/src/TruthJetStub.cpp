// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This is a stub version when the FastJet plugin is not available
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace Acts::Examples::Python {
void exportTruthJet(py::module& m) {
  // Empty implementation when FastJet is not available
}
}  // namespace Acts::Examples::Python
