// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Test-only pybind11 module (NOT part of the production `acts` package). It
// provides a Core-only way to disown a smart_holder container the same way the
// whiteboard does (by taking it as a std::unique_ptr<T>), so the fail-loud
// behaviour of the EventData proxy tethers can be tested without depending on
// acts.examples. pybind11's type registry is process-wide, so these functions
// resolve the SpacePointContainer2 / SeedContainer2 types registered by
// `import acts`.

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(_acts_core_test_bindings, m) {
  m.doc() = "Test-only helpers to disown ACTS EventData containers";

  // Taking the argument by std::unique_ptr<T> tells pybind11's smart_holder to
  // disown the Python wrapper (and frees the object when the unique_ptr goes
  // out of scope) -- exactly what WhiteBoardRegistry::fromPython does.
  m.def("consume_spacepoints",
        [](std::unique_ptr<Acts::SpacePointContainer2>) {});
  m.def("consume_seeds", [](std::unique_ptr<Acts::SeedContainer2>) {});
}
