// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "ActsPlugins/Detray/DetrayGeometryConverter.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
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

  using DetrayMetaDataODD = detray::odd_metadata<detray::array<float>>;
  using DetrayDetectorODD = detray::detector<DetrayMetaDataODD>;

  py::class_<DetrayDetectorODD::name_map>(detray, "DetrayDetectorODDNameMap");

  // Register the raw detector type so shared_ptr<DetrayDetectorODD> can be
  // passed to functions in other pybind11 modules (e.g.
  // StraightLinePropagatorODD).
  py::class_<DetrayDetectorODD, std::shared_ptr<DetrayDetectorODD>>(
      detray, "DetrayDetectorODD")
      .def("volumes", [](DetrayDetectorODD& self) { return self.volumes(); })
      .def("surfaces", [](DetrayDetectorODD& self) { return self.surfaces(); })
      .def("checkConsistency",
           [](DetrayDetectorODD& self) {
             detray::detail::check_consistency(self);
           })
      .def("writeToJson",
           [](DetrayDetectorODD& self, const DetrayDetectorODD::name_map& names,
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

  // ── Payload converter ──────────────────────────────────────────────────
  auto payloadConverter = py::class_<DetrayPayloadConverter,
                                     std::shared_ptr<DetrayPayloadConverter>>(
      detray, "DetrayPayloadConverter");

  auto payloadConfig =
      py::class_<DetrayPayloadConverter::Config>(payloadConverter, "Config")
          .def(py::init<>())
          .def_readwrite("sensitiveStrategy",
                         &DetrayPayloadConverter::Config::sensitiveStrategy)
          .def_property(
              "beampipeVolume",
              [](const DetrayPayloadConverter::Config& cfg) {
                return cfg.beampipeVolume;
              },
              [](DetrayPayloadConverter::Config& cfg,
                 const Acts::TrackingVolume* volume) {
                cfg.beampipeVolume = volume;
              },
              py::return_value_policy::reference);

  py::enum_<DetrayPayloadConverter::Config::SensitiveStrategy>(
      payloadConfig, "SensitiveStrategy")
      .value("Identifier",
             DetrayPayloadConverter::Config::SensitiveStrategy::Identifier)
      .value(
          "DetectorElement",
          DetrayPayloadConverter::Config::SensitiveStrategy::DetectorElement);

  payloadConverter.def(
      py::init([](const DetrayPayloadConverter::Config& config,
                  Acts::Logging::Level level) {
        return std::make_shared<DetrayPayloadConverter>(
            config, Acts::getDefaultLogger("DetrayPayloadConverter", level));
      }),
      "config"_a, "level"_a = Acts::Logging::INFO);

  // ── Geometry converter ─────────────────────────────────────────────────
  auto geometryConverter = py::class_<DetrayGeometryConverter,
                                      std::shared_ptr<DetrayGeometryConverter>>(
      detray, "DetrayGeometryConverter");

  py::class_<DetrayGeometryConverter::Config>(geometryConverter, "Config")
      .def(py::init<>())
      .def_readwrite("payloadConverter",
                     &DetrayGeometryConverter::Config::payloadConverter)
      .def_readwrite("convertMaterial",
                     &DetrayGeometryConverter::Config::convertMaterial)
      .def_readwrite("convertSurfaceGrids",
                     &DetrayGeometryConverter::Config::convertSurfaceGrids);

  using DetrayGeometryODD =
      DetrayGeometryConverter::DetrayGeometry<DetrayMetaDataODD>;
  py::class_<DetrayGeometryODD>(geometryConverter, "DetrayGeometry")
      .def_readonly("detector", &DetrayGeometryODD::detector)
      .def_readonly("names", &DetrayGeometryODD::names);

  geometryConverter
      .def(py::init([](DetrayGeometryConverter::Config config,
                       Acts::Logging::Level level) {
             return std::make_shared<DetrayGeometryConverter>(
                 std::move(config),
                 Acts::getDefaultLogger("DetrayGeometryConverter", level));
           }),
           "config"_a, "level"_a = Acts::Logging::INFO)
      .def(
          "convert",
          [](const DetrayGeometryConverter& self, vecmem::memory_resource& mr,
             const Acts::GeometryContext& gctx,
             std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
             const std::string& detectorName) {
            return self.convert<DetrayMetaDataODD>(
                mr, gctx, std::move(trackingGeometry), detectorName);
          },
          "mr"_a, "gctx"_a, "trackingGeometry"_a, "detectorName"_a = "");
}
