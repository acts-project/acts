// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "ActsPlugins/Detray/DetrayDetectorIO.hpp"
#include "ActsPlugins/Detray/DetrayGeometryConverter.hpp"
#include "ActsPlugins/Detray/DetrayMetadata.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <detray/core/detector.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsDetray, detray) {
  using namespace ActsPlugins;

  // The detray volume/surface name map is the same type for every metadata, so
  // it is registered exactly once.
  py::class_<detray::name_map>(detray, "DetrayNameMap");

  // Metadata markers for the closed set. These empty classes are passed (as the
  // class object) to DetrayGeometryConverter.convert to select the metadata,
  // mirroring the C++ template argument convert<Metadata>.
  auto oddMetadata = py::class_<DetrayMetadata::Odd>(detray, "OddMetadata");
  auto defaultMetadata =
      py::class_<DetrayMetadata::Default>(detray, "DefaultMetadata");

  // Read a pre-built detray detector from JSON file(s). The metadata is
  // selected by passing one of the metadata markers, mirroring convert.
  detray.def(
      "read",
      [oddType = py::object(oddMetadata),
       defaultType = py::object(defaultMetadata)](
          py::object metadata, vecmem::memory_resource& mr,
          const std::vector<std::string>& files) -> py::object {
        if (metadata.is(oddType)) {
          return py::cast(
              detail::readDetrayDetector<DetrayMetadata::Odd>(mr, files));
        }
        if (metadata.is(defaultType)) {
          return py::cast(
              detail::readDetrayDetector<DetrayMetadata::Default>(mr, files));
        }
        throw std::invalid_argument(
            "detray.read: unsupported metadata; pass one of the metadata "
            "markers (e.g. acts.detray.OddMetadata)");
      },
      "metadata"_a, "mr"_a, "files"_a);

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

  // Register the per-metadata detector and conversion-result types.
  auto registerMetadata = [&](auto metadataTag, const std::string& suffix) {
    using Metadata = typename decltype(metadataTag)::type;
    using Detector = detray::detector<Metadata>;
    using Geometry = DetrayGeometryConverter::DetrayGeometry<Metadata>;

    py::class_<Detector, std::shared_ptr<Detector>>(
        detray, ("DetrayDetector" + suffix).c_str())
        .def("volumes", [](Detector& self) { return self.volumes(); })
        .def("surfaces", [](Detector& self) { return self.surfaces(); })
        .def("checkConsistency",
             [](Detector& self) { detail::checkDetrayConsistency(self); })
        .def("writeToJson",
             [](Detector& self, const detray::name_map& names,
                const std::string& fname) {
               detail::writeDetrayJson(self, names, fname);
             });

    py::class_<Geometry>(geometryConverter, ("DetrayGeometry" + suffix).c_str())
        .def_readonly("detector", &Geometry::detector)
        .def_readonly("names", &Geometry::names);
  };

  registerMetadata(std::type_identity<DetrayMetadata::Odd>{}, "ODD");
  registerMetadata(std::type_identity<DetrayMetadata::Default>{}, "Default");

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
          [oddType = py::object(oddMetadata),
           defaultType = py::object(defaultMetadata)](
              const DetrayGeometryConverter& self, py::object metadata,
              vecmem::memory_resource& mr, const Acts::GeometryContext& gctx,
              std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
              const std::string& detectorName) -> py::object {
            if (metadata.is(oddType)) {
              return py::cast(self.convert<DetrayMetadata::Odd>(
                  mr, gctx, std::move(trackingGeometry), detectorName));
            }
            if (metadata.is(defaultType)) {
              return py::cast(self.convert<DetrayMetadata::Default>(
                  mr, gctx, std::move(trackingGeometry), detectorName));
            }
            throw std::invalid_argument(
                "DetrayGeometryConverter.convert: unsupported metadata; pass "
                "one of the metadata markers (e.g. acts.detray.OddMetadata)");
          },
          "metadata"_a, "mr"_a, "gctx"_a, "trackingGeometry"_a,
          "detectorName"_a = "");
}
