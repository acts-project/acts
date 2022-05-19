// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimHitReader.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace Acts::Python {
void addEDM4hep(Context& ctx) {
  auto mex = ctx.get("examples");

  {
    using Reader = ActsExamples::EDM4hepSimHitReader;
    using Config = Reader::Config;
    auto reader =
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
            mex, "EDM4hepSimHitReader")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(reader, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputPath);
    ACTS_PYTHON_MEMBER(outputSimHits);
    ACTS_PYTHON_MEMBER(dd4hepGeometryService);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python
