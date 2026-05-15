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

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelDetector.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelMuonMockupBuilder.hpp"
#include "ActsExamples/ITkModuleSplitting/ITkModuleSplitting.hpp"
#include "ActsExamples/MuonSpectrometerMockupDetector/GeoMuonMockupExperiment.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <string>

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoVPhysVol.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsPlugins;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsGeoModel, gm) {
  {
    auto f =
        py::class_<GeoModelDetector, Detector,
                   std::shared_ptr<GeoModelDetector>>(gm, "GeoModelDetector")
            .def(py::init<const GeoModelDetector::Config&>())
            .def("buildTrackingGeometry",
                 &GeoModelDetector::buildTrackingGeometry);

    auto c =
        py::class_<GeoModelDetector::Config>(f, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, geoModelTree, logLevel, path);

    // patch the constructor
    patchKwargsConstructor(c);
  }
  {
    // GeomodelMuonMockupBuilder
    py::class_<ITrackingGeometryBuilder,
               std::shared_ptr<ITrackingGeometryBuilder>>(
        gm, "ITrackingGeometryBuilder");
    auto gmMuonBuilder =
        py::class_<GeoModelMuonMockupBuilder, ITrackingGeometryBuilder,
                   std::shared_ptr<GeoModelMuonMockupBuilder>>(
            gm, "GeoModelMuonMockupBuilder")
            .def(py::init([](const GeoModelMuonMockupBuilder::Config& config,
                             const std::string& name, Logging::Level level) {
              return std::make_shared<GeoModelMuonMockupBuilder>(
                  config, getDefaultLogger(name, level));
            }))
            .def("trackingGeometry",
                 &GeoModelMuonMockupBuilder::trackingGeometry);
    auto gmMuonConfig =
        py::class_<GeoModelMuonMockupBuilder::Config>(gmMuonBuilder, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(gmMuonConfig, volumeBoxFPVs, stationNames,
                       volumeBoundFactory);
  }

  {
    // GeoMuonMockupExperiment
    auto f =
        py::class_<GeoMuonMockupExperiment,
                   std::shared_ptr<GeoMuonMockupExperiment>>(
            gm, "GeoMuonMockupExperiment")
            .def(py::init([](const GeoMuonMockupExperiment::Config& config,
                             const std::string& name, Logging::Level level) {
              return std::make_shared<GeoMuonMockupExperiment>(
                  config, getDefaultLogger(name, level));
            }))
            .def("constructMS", &GeoMuonMockupExperiment::constructMS);
    auto c = py::class_<GeoMuonMockupExperiment::Config>(f, "Config")
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
                       buildEndcaps, nInnerMultiplets, nInnerGasGapsPerMl);
  }

  /// Module splitting
  {
    gm.def(
        "splitBarrelModule",
        [](const GeometryContext& gctx,
           std::shared_ptr<const GeoModelDetectorElement> detElement,
           unsigned nSegments, Logging::Level logLevel) {
          auto logger = getDefaultLogger("ITkSlitBarrel", logLevel);
          auto name = detElement->databaseEntryName();

          auto factory = [&](const auto& trafo, const auto& bounds) {
            return GeoModelDetectorElement::createDetectorElement<
                PlaneSurface, RectangleBounds>(detElement->physicalVolume(),
                                               bounds, trafo,
                                               detElement->thickness());
          };

          return ITk::splitBarrelModule(gctx, detElement, nSegments, factory,
                                        name, *logger);
        },
        "gxtx"_a, "detElement"_a, "nSegments"_a, "logLevel"_a = Logging::INFO);

    gm.def(
        "splitDiscModule",
        [](const GeometryContext& gctx,
           std::shared_ptr<const GeoModelDetectorElement> detElement,
           const std::vector<std::pair<double, double>>& patterns,
           Logging::Level logLevel) {
          auto logger = getDefaultLogger("ITkSlitBarrel", logLevel);
          auto name = detElement->databaseEntryName();

          auto factory = [&](const auto& trafo, const auto& bounds) {
            return GeoModelDetectorElement::createDetectorElement<
                DiscSurface, AnnulusBounds>(detElement->physicalVolume(),
                                            bounds, trafo,
                                            detElement->thickness());
          };

          return ITk::splitDiscModule(gctx, detElement, patterns, factory, name,
                                      *logger);
        },
        "gxtx"_a, "detElement"_a, "splitRanges"_a,
        "logLevel"_a = Logging::INFO);
  }
}
