// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Detray/DetrayGeometryConverter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <string>

#include <detray/builders/detector_builder.hpp>
#include <detray/detectors/default_metadata.hpp>
#include <detray/detectors/odd_metadata.hpp>
#include <pybind11/pybind11.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsDetray, detray) {
  using namespace ActsPlugins;

  {
    using DetrayMetaDataODD = detray::odd_metadata<detray::array<float>>;

    using DetrayDetectorODD = detray::detector<DetrayMetaDataODD>;

    py::class_<DetrayDetectorODD, std::shared_ptr<DetrayDetectorODD>>(
        detray, "DetrayDetectorODD")
        .def("volumes", &DetrayDetectorODD::volumes)
        .def("surfaces", &DetrayDetectorODD::surfaces)
        .def("check_consistency", [](DetrayDetectorODD& self) {
          detray::detail::check_consistency(self);
        });

    detray.def(
        "convertODD",
        [](vecmem::memory_resource& mr, const Acts::GeometryContext& gctx,
           const Acts::TrackingGeometry& trackingGeometry,
           const std::string& beampipeVolumeName,
           Acts::Logging::Level logLevel = Acts::Logging::INFO) {
          return DetrayGeometryConverter::toDetray<DetrayMetaDataODD>(
              mr, gctx, trackingGeometry, beampipeVolumeName, logLevel);
        },
        "mr"_a, "gctx"_a, "trackingGeometry"_a, "beampipeVolumeName"_a,
        "logLevel"_a = Acts::Logging::INFO);
  }
}
