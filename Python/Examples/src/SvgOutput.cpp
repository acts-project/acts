// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Io/Svg/SvgPointWriter.hpp"
#include "ActsExamples/Io/Svg/SvgTrackingGeometryWriter.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <algorithm>
#include <memory>
#include <ranges>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace ActsPython {
/// This adds the SVG output related components to the examples module
/// @param mex the examples module
void addSvgOutput(py::module_& mex) {
  // Components from the ActsExamples - part of acts.examples
  {
    using Writer = SvgTrackingGeometryWriter;
    auto w = py::class_<Writer, std::shared_ptr<Writer>>(
                 mex, "SvgTrackingGeometryWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", py::overload_cast<const AlgorithmContext&,
                                                 const Acts::TrackingGeometry&>(
                                   &Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, outputDir, converterOptions);
  }
  {
    using Writer = SvgPointWriter<SimSpacePoint>;
    auto w =
        py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
            mex, "SvgSimSpacePointWriter")
            .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write",
                 py::overload_cast<const AlgorithmContext&>(&Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, writerName, trackingGeometry, inputCollection,
                       infoBoxTitle, outputDir);
  }

  {
    using Writer = SvgPointWriter<SimHit, AccessorPositionXYZ>;
    auto w =
        py::class_<Writer, IWriter, std::shared_ptr<Writer>>(mex,
                                                             "SvgSimHitWriter")
            .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write",
                 py::overload_cast<const AlgorithmContext&>(&Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, writerName, trackingGeometry, inputCollection,
                       infoBoxTitle, outputDir);
  }
}
}  // namespace ActsPython
