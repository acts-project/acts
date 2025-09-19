// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

#include <array>
#include <exception>
#include <memory>
#include <string>
#include <unordered_map>

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Acts;

namespace ActsPython {

void addContext(Context& ctx) {
  auto& m = ctx.get("main");

  py::class_<Acts::MagneticFieldContext>(m, "MagneticFieldContext")
      .def(py::init<>());
}

}  // namespace ActsPython
