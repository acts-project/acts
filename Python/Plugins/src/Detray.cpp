// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "ActsPlugins/Detray/DetrayGeometryConverter.hpp"
#include "ActsPlugins/Detray/DetrayJsonConverter.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
#include "ActsPlugins/Detray/DetraySurfaceGridPolicy.hpp"
#include "ActsPlugins/Detray/DetraySurfaceIdMap.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <optional>
#include <string>

#include <detray/io/frontend/detector_reader_config.hpp>
#include <detray/io/frontend/detector_writer.hpp>
#include <detray/io/frontend/detector_writer_config.hpp>
#include <detray/io/json/detector_reader.hpp>
#include <detray/python/detector.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsDetray, detray) {
  using namespace ActsPlugins;

  // Detray owns the canonical detector, allocator and name-map registrations.
  auto core = py::module_::import("detray.core");
  detray.attr("Detector") = core.attr("Detector");

  // ── Surface grid policy ────────────────────────────────────────────────
  {
    py::class_<DetraySurfaceGridConfig>(detray, "DetraySurfaceGridConfig")
        .def(py::init<>())
        .def_readwrite("vertexZHalfRange",
                       &DetraySurfaceGridConfig::vertexZHalfRange)
        .def_readwrite("minPt", &DetraySurfaceGridConfig::minPt)
        .def_readwrite("bField", &DetraySurfaceGridConfig::bField)
        .def_readwrite("nVertexZ", &DetraySurfaceGridConfig::nVertexZ)
        .def_readwrite("nCurvature", &DetraySurfaceGridConfig::nCurvature)
        .def_readwrite("maxSampleSpacing",
                       &DetraySurfaceGridConfig::maxSampleSpacing)
        .def_readwrite("numberOfBinsFactor",
                       &DetraySurfaceGridConfig::numberOfBinsFactor);

    py::class_<DetraySurfaceGridPolicy>(detray, "DetraySurfaceGridPolicy")
        .def_property_readonly("config", &DetraySurfaceGridPolicy::config);

    detray.def(
        "addDetraySurfaceGridPolicy",
        [](Acts::NavigationPolicyFactory& factory,
           const DetraySurfaceGridConfig& config) {
          return std::move(factory).template add<DetraySurfaceGridPolicy>(
              config);
        },
        "factory"_a, "config"_a = DetraySurfaceGridConfig{},
        "Add a DetraySurfaceGridPolicy carrying @p config to the factory");

    detray.def("addDetrayNavigationPolicies", &addDetrayNavigationPolicies,
               "config"_a,
               "Register the Detray navigation policies with a "
               "TrackingGeometryJsonConverter config");
  }

  // ── Payload converter ──────────────────────────────────────────────────
  auto payloadConverter = py::class_<DetrayPayloadConverter,
                                     std::shared_ptr<DetrayPayloadConverter>>(
      detray, "DetrayPayloadConverter");

  auto payloadConfig =
      py::class_<DetrayPayloadConverter::Config>(payloadConverter, "Config")
          .def(py::init<>())
          .def_readwrite("sensitiveStrategy",
                         &DetrayPayloadConverter::Config::sensitiveStrategy)
          .def_readwrite("includeWorld",
                         &DetrayPayloadConverter::Config::includeWorld)
          .def_readwrite(
              "portalSegmentTolerance",
              &DetrayPayloadConverter::Config::portalSegmentTolerance)
          .def_readwrite("convertMaterial",
                         &DetrayPayloadConverter::Config::convertMaterial)
          .def_readwrite("convertSurfaceGrids",
                         &DetrayPayloadConverter::Config::convertSurfaceGrids)
          .def_property(
              "detectorName",
              [](const DetrayPayloadConverter::Config& cfg) {
                return cfg.detectorName;
              },
              [](DetrayPayloadConverter::Config& cfg, const std::string& name) {
                cfg.detectorName = name;
              },
              py::return_value_policy::reference)
          .def_property(
              "beampipeVolume",
              [](const DetrayPayloadConverter::Config& cfg) {
                return cfg.beampipeVolume;
              },
              [](DetrayPayloadConverter::Config& cfg,
                 const Acts::TrackingVolume* volume) {
                cfg.beampipeVolume = volume;
              },
              py::return_value_policy::reference)
          .def_property(
              "replacementWorld",
              [](const DetrayPayloadConverter::Config& cfg) {
                return cfg.replacementWorld;
              },
              [](DetrayPayloadConverter::Config& cfg,
                 const Acts::TrackingVolume* volume) {
                cfg.replacementWorld = volume;
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
      .def_readwrite("deduplicateMaterial",
                     &DetrayGeometryConverter::Config::deduplicateMaterial);

  using DetrayGeometry = DetrayGeometryConverter::RegisteredDetrayGeometry;
  py::class_<DetraySurfaceIdMap, std::shared_ptr<DetraySurfaceIdMap>>(
      detray, "DetraySurfaceIdMap")
      // detray identifiers are exposed as their integer value
      .def(
          "toDetray",
          [](const DetraySurfaceIdMap& self, Acts::GeometryIdentifier id)
              -> std::optional<detray::geometry::identifier::value_t> {
            if (auto detrayId = self.toDetray(id)) {
              return detrayId->value();
            }
            return std::nullopt;
          },
          "id"_a)
      .def(
          "fromDetray",
          [](const DetraySurfaceIdMap& self,
             detray::geometry::identifier::value_t detrayId) {
            return self.fromDetray(detray::geometry::identifier{detrayId});
          },
          "detrayId"_a)
      .def("__len__", &DetraySurfaceIdMap::size);

  py::class_<DetrayGeometry>(geometryConverter, "DetrayGeometry")
      .def_property_readonly("detector",
                             [](py::object owner) {
                               auto& result = owner.cast<DetrayGeometry&>();
                               py::object detector = py::cast(result.detector);
                               // The result retains the allocator and ACTS
                               // geometry. Preserve that ownership chain when
                               // callers extract just the detector.
                               py::detail::keep_alive_impl(detector, owner);
                               return detector;
                             })
      .def_readonly("names", &DetrayGeometry::names)
      .def_readonly("surfaceIds", &DetrayGeometry::surfaceIds);

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
             const std::string& metadataId) {
            return self.convertRegistered<detray::python::algebra>(
                mr, gctx, trackingGeometry, metadataId);
          },
          "mr"_a, "gctx"_a, "trackingGeometry"_a, "metadataId"_a = "odd",
          py::keep_alive<0, 2>(), py::keep_alive<0, 4>());
}
