// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Must be on top to avoid some conflict between forward declare and typedef
// Needed until https://gitlab.cern.ch/GeoModelDev/GeoModel/-/merge_requests/351
// is deployed
// clang-format off
#include <GeoModelRead/ReadGeoModel.h>
// clang-format on

#include "Acts/Detector/CylindricalContainerBuilder.hpp"
#include "Acts/Geometry/ITrackingGeometryBuilder.hpp"
#include "Acts/Plugins/GeoModel/GeoModelBlueprintCreater.hpp"
#include "Acts/Plugins/GeoModel/GeoModelConverters.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElementITk.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Plugins/GeoModel/IGeoShapeConverter.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelDetector.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelMuonMockupBuilder.hpp"
#include "ActsExamples/ITkModuleSplitting/ITkModuleSplitting.hpp"
#include "ActsExamples/MuonSpectrometerMockupDetector/GeoMuonMockupExperiment.hpp"
#include "ActsPython/Utilities/Context.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include "ActsPython/Utilities/Patchers.hpp"

#include <string>

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoVPhysVol.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace ActsPython {

void addGeoModelDetector(Context& ctx) {
  auto& mex = ctx.get("mex");

  auto gm = mex.def_submodule("geomodel");

  {
    auto f = py::class_<ActsExamples::GeoModelDetector, ActsExamples::Detector,
                        std::shared_ptr<ActsExamples::GeoModelDetector>>(
                 gm, "GeoModelDetector")
                 .def(py::init<const ActsExamples::GeoModelDetector::Config&>())
                 .def("buildTrackingGeometry",
                      &ActsExamples::GeoModelDetector::buildTrackingGeometry);

    auto c = py::class_<ActsExamples::GeoModelDetector::Config>(f, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, geoModelTree, logLevel, path);

    // patch the constructor
    patchKwargsConstructor(c);
  }
  {
    // GeomodelMuonMockupBuilder
    py::class_<Acts::ITrackingGeometryBuilder,
               std::shared_ptr<Acts::ITrackingGeometryBuilder>>(
        gm, "ITrackingGeometryBuilder");
    auto gmMuonBuilder =
        py::class_<ActsExamples::GeoModelMuonMockupBuilder,
                   Acts::ITrackingGeometryBuilder,
                   std::shared_ptr<ActsExamples::GeoModelMuonMockupBuilder>>(
            gm, "GeoModelMuonMockupBuilder")
            .def(py::init([](const ActsExamples::GeoModelMuonMockupBuilder::
                                 Config& config,
                             const std::string& name,
                             Acts::Logging::Level level) {
              return std::make_shared<ActsExamples::GeoModelMuonMockupBuilder>(
                  config, Acts::getDefaultLogger(name, level));
            }))
            .def("trackingGeometry",
                 &ActsExamples::GeoModelMuonMockupBuilder::trackingGeometry);
    auto gmMuonConfig =
        py::class_<ActsExamples::GeoModelMuonMockupBuilder::Config>(
            gmMuonBuilder, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(gmMuonConfig, volumeBoxFPVs, stationNames,
                       volumeBoundFactory);
  }
  {
    // GeoMuonMockupExperiment
    auto f =
        py::class_<ActsExamples::GeoMuonMockupExperiment,
                   std::shared_ptr<ActsExamples::GeoMuonMockupExperiment>>(
            gm, "GeoMuonMockupExperiment")
            .def(py::init(
                [](const ActsExamples::GeoMuonMockupExperiment::Config& config,
                   const std::string& name, Acts::Logging::Level level) {
                  return std::make_shared<
                      ActsExamples::GeoMuonMockupExperiment>(
                      config, Acts::getDefaultLogger(name, level));
                }))
            .def("constructMS",
                 &ActsExamples::GeoMuonMockupExperiment::constructMS);
    auto c =
        py::class_<ActsExamples::GeoMuonMockupExperiment::Config>(f, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(c,
                       /// General properties
                       dumpTree, dbName,
                       /// Mdt properties
                       innerTubeRadius, tubeWallThickness, nTubeLayers, nTubes,
                       mdtFoamThickness, multiLayerSeparation,
                       /// Rpc properties
                       nRpcGasGaps, nRpcAlongZ, nRpcAlongPhi,
                       /// Station properties
                       barrelRadii, nSectors, nEtaStations, stationDistInZ,
                       stationDistInR, endCapWheelLowR, bigWheelDistZ,
                       buildEndcaps);
  }

  gm.def(
      "splitBarrelModule",
      [](const Acts::GeometryContext& gctx,
         std::shared_ptr<const Acts::GeoModelDetectorElement> detElement,
         unsigned nSegments, Acts::Logging::Level logLevel) {
        auto logger = Acts::getDefaultLogger("ITkSlitBarrel", logLevel);
        auto name = detElement->databaseEntryName();

        auto factory = [&](const auto& trafo, const auto& bounds) {
          return Acts::GeoModelDetectorElement::createDetectorElement<
              Acts::PlaneSurface, Acts::RectangleBounds>(
              detElement->physicalVolume(), bounds, trafo,
              detElement->thickness());
        };

        return ActsExamples::ITk::splitBarrelModule(gctx, detElement, nSegments,
                                                    factory, name, *logger);
      },
      "gxtx"_a, "detElement"_a, "nSegments"_a,
      "logLevel"_a = Acts::Logging::INFO);

  gm.def(
      "splitDiscModule",
      [](const Acts::GeometryContext& gctx,
         std::shared_ptr<const Acts::GeoModelDetectorElement> detElement,
         const std::vector<std::pair<double, double>>& patterns,
         Acts::Logging::Level logLevel) {
        auto logger = Acts::getDefaultLogger("ITkSlitBarrel", logLevel);
        auto name = detElement->databaseEntryName();

        auto factory = [&](const auto& trafo, const auto& bounds) {
          return Acts::GeoModelDetectorElement::createDetectorElement<
              Acts::DiscSurface, Acts::AnnulusBounds>(
              detElement->physicalVolume(), bounds, trafo,
              detElement->thickness());
        };

        return ActsExamples::ITk::splitDiscModule(gctx, detElement, patterns,
                                                  factory, name, *logger);
      },
      "gxtx"_a, "detElement"_a, "splitRanges"_a,
      "logLevel"_a = Acts::Logging::INFO);
}

}  // namespace ActsPython
