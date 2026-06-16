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
#include <detray/io/frontend/detector_reader.hpp>
#include <detray/io/frontend/detector_reader_config.hpp>
#include <detray/io/frontend/detector_writer.hpp>
#include <detray/io/frontend/detector_writer_config.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsDetray, detray) {
  using namespace ActsPlugins;

  {
    using DetrayMetaDataODD = detray::odd_metadata<detray::array<float>>;
    using DetrayDetectorODD = detray::detector<DetrayMetaDataODD>;

    py::class_<DetrayDetectorODD::name_map>(detray, "DetrayDetectorODDNameMap");

    // Register the raw detector type so shared_ptr<DetrayDetectorODD> can be
    // passed to functions in other pybind11 modules (e.g.
    // StraightLinePropagatorODD).
    py::class_<DetrayDetectorODD, std::shared_ptr<DetrayDetectorODD>>(
        detray, "DetrayDetectorODD")
        .def("volumes", [](DetrayDetectorODD& self) { return self.volumes(); })
        .def("surfaces",
             [](DetrayDetectorODD& self) { return self.surfaces(); })
        .def("checkConsistency",
             [](DetrayDetectorODD& self) {
               detray::detail::check_consistency(self);
             })
        .def("writeToJson", [](DetrayDetectorODD& self,
                               const DetrayDetectorODD::name_map& names,
                               const std::string& fname) {
          auto cfg = detray::io::detector_writer_config{}
                         .format(detray::io::format::json)
                         .path(fname)
                         .replace_files(true);
          detray::io::write_detector(self, names, cfg);
        });

    detray.def(
        "readODD",
        [](vecmem::memory_resource& mr, const std::vector<std::string>& files) {
          auto cfg = detray::io::detector_reader_config{}.do_check(false);
          for (const auto& f : files) {
            cfg.add_file(f);
          }
          return detray::io::read_detector<DetrayDetectorODD>(mr, cfg);
        },
        "mr"_a, "files"_a);

    detray.def(
        "convertODD",
        [](vecmem::memory_resource& mr, const Acts::GeometryContext& gctx,
           const Acts::TrackingGeometry& trackingGeometry,
           const std::string& beampipeVolumeName,
           const std::string& detectorName = "",
           Acts::Logging::Level logLevel = Acts::Logging::INFO,
           bool convertMaterial = true, bool convertSurfaceGrids = true) {
          auto [det, names] =
              DetrayGeometryConverter::toDetray<DetrayMetaDataODD>(
                  mr, gctx, trackingGeometry, beampipeVolumeName, detectorName,
                  logLevel, convertMaterial, convertSurfaceGrids);
          return std::make_pair(std::move(det), std::move(names));
        },
        "mr"_a, "gctx"_a, "trackingGeometry"_a, "beampipeVolumeName"_a,
        "detectorName"_a = "", "logLevel"_a = Acts::Logging::INFO,
        "convertMaterial"_a = true, "convertSurfaceGrids"_a = true);
  }
}
