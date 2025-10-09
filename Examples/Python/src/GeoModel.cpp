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
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelDetector.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelMuonMockupBuilder.hpp"
#include "ActsExamples/ITkModuleSplitting/ITkModuleSplitting.hpp"
#include "ActsExamples/MuonSpectrometerMockupDetector/GeoMuonMockupExperiment.hpp"
#include "ActsPlugins/GeoModel/GeoModelBlueprintCreater.hpp"
#include "ActsPlugins/GeoModel/GeoModelConverters.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorElement.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorElementITk.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
#include "ActsPlugins/GeoModel/GeoModelReader.hpp"
#include "ActsPlugins/GeoModel/GeoModelTree.hpp"
#include "ActsPlugins/GeoModel/IGeoShapeConverter.hpp"
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

namespace ActsPython {

void addGeoModel(Context& ctx) {
  auto m = ctx.get("main");

  auto gm = m.def_submodule("geomodel");

  py::class_<GeoModelTree::FpvConstLink>(gm, "GeoModelTree::FpvConstLink")
      .def(py::init<>())
      .def("get", &GeoModelTree::FpvConstLink::get,
           py::return_value_policy::reference);

  py::class_<GeoModelTree>(gm, "GeoModelTree").def(py::init<>());

  gm.def("readFromDb", &GeoModelReader::readFromDb);

  py::class_<GeoModelDetectorElement, std::shared_ptr<GeoModelDetectorElement>>(
      gm, "GeoModelDetectorElement")
      .def("logVolName", &GeoModelDetectorElement::logVolName)
      .def("databaseEntryName", &GeoModelDetectorElement::databaseEntryName)
      .def("surface", [](GeoModelDetectorElement self) {
        return self.surface().getSharedPtr();
      });

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

  // Shape converters
  {
    py::class_<IGeoShapeConverter, std::shared_ptr<IGeoShapeConverter>>(
        gm, "IGeoShapeConverter");

    py::class_<GeoBoxConverter, IGeoShapeConverter,
               std::shared_ptr<GeoBoxConverter>>(gm, "GeoBoxConverter")
        .def(py::init<>())
        .def("toSensitiveSurface", &GeoBoxConverter::toSensitiveSurface)
        .def("toPassiveSurface", &GeoBoxConverter::toPassiveSurface);

    py::class_<GeoTrdConverter, IGeoShapeConverter,
               std::shared_ptr<GeoTrdConverter>>(gm, "GeoTrdConverter")
        .def(py::init<>())
        .def("toSensitiveSurface", &GeoTrdConverter::toSensitiveSurface)
        .def("toPassiveSurface", &GeoTrdConverter::toPassiveSurface);

    py::class_<GeoTubeConverter, IGeoShapeConverter,
               std::shared_ptr<GeoTubeConverter>>(gm, "GeoTubeConverter")
        .def(py::init<>())
        .def("toSensitiveSurface", &GeoTubeConverter::toSensitiveSurface)
        .def("toPassiveSurface", &GeoTubeConverter::toPassiveSurface);

    py::class_<GeoUnionDoubleTrdConverter, IGeoShapeConverter,
               std::shared_ptr<GeoUnionDoubleTrdConverter>>(
        gm, "GeoUnionDoubleTrdConverter")
        .def(py::init<>())
        .def("toSensitiveSurface",
             &GeoUnionDoubleTrdConverter::toSensitiveSurface)
        .def("toPassiveSurface", &GeoUnionDoubleTrdConverter::toPassiveSurface);

    py::class_<GeoIntersectionAnnulusConverter, IGeoShapeConverter,
               std::shared_ptr<GeoIntersectionAnnulusConverter>>(
        gm, "GeoIntersectionAnnulusConverter")
        .def(py::init<>())
        .def("toSensitiveSurface",
             &GeoIntersectionAnnulusConverter::toSensitiveSurface)
        .def("toPassiveSurface",
             &GeoIntersectionAnnulusConverter::toPassiveSurface);

    py::class_<GeoShiftConverter, IGeoShapeConverter,
               std::shared_ptr<GeoShiftConverter>>(gm, "GeoShiftConverter")
        .def(py::init<>())
        .def("toSensitiveSurface", &GeoShiftConverter::toSensitiveSurface)
        .def("toPassiveSurface", &GeoShiftConverter::toPassiveSurface);
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
  // Volume factory
  {
    auto a =
        py::class_<GeoModelDetectorObjectFactory,
                   std::shared_ptr<GeoModelDetectorObjectFactory>>(
            gm, "GeoModelDetectorObjectFactory")
            .def(py::init([](const GeoModelDetectorObjectFactory::Config& cfg,
                             Logging::Level level) {
              return std::make_shared<GeoModelDetectorObjectFactory>(
                  cfg,
                  getDefaultLogger("GeoModelDetectorObjectFactory", level));
            }))
            .def("construct", &GeoModelDetectorObjectFactory::construct);

    py::class_<GeoModelDetectorObjectFactory::Config>(a, "Config")
        .def(py::init<>())
        .def_readwrite(
            "convertSubVolumes",
            &GeoModelDetectorObjectFactory::Config::convertSubVolumes)
        .def_readwrite("nameList",
                       &GeoModelDetectorObjectFactory::Config::nameList)
        .def_readwrite("convertBox",
                       &GeoModelDetectorObjectFactory::Config::convertBox)
        .def_readwrite("materialList",
                       &GeoModelDetectorObjectFactory::Config::materialList);

    auto convVol = py::class_<GeoModelDetectorObjectFactory::ConvertedGeoVol>(
        a, "ConvertedGeoVol");

    ACTS_PYTHON_STRUCT(convVol, volume, gen2Volume, fullPhysVol, name,
                       surfaces);
    py::class_<GeoModelDetectorObjectFactory::Cache>(a, "Cache")
        .def(py::init<>())
        .def_readwrite("sensitiveSurfaces",
                       &GeoModelDetectorObjectFactory::Cache::sensitiveSurfaces)
        .def_readwrite("boundingBoxes",
                       &GeoModelDetectorObjectFactory::Cache::volumeBoxFPVs);

    py::class_<GeoModelDetectorObjectFactory::Options>(a, "Options")
        .def(py::init<>())
        .def_readwrite("queries",
                       &GeoModelDetectorObjectFactory::Options::queries);
  }

  {
    py::class_<GeoModelBlueprintCreater::Blueprint,
               std::shared_ptr<GeoModelBlueprintCreater::Blueprint>>(
        gm, "Blueprint")
        .def("convertToBuilder", [](GeoModelBlueprintCreater::Blueprint& self,
                                    Logging::Level level) {
          // It's a container builder
          return std::make_shared<Experimental::CylindricalContainerBuilder>(
              self.node(), level);
        });

    auto bpc =
        py::class_<GeoModelBlueprintCreater,
                   std::shared_ptr<GeoModelBlueprintCreater>>(
            gm, "GeoModelBlueprintCreater")
            .def(py::init([](const GeoModelBlueprintCreater::Config& cfg,
                             Logging::Level level) {
              return std::make_shared<GeoModelBlueprintCreater>(
                  cfg, getDefaultLogger("GeoModelBlueprintCreater", level));
            }))
            .def("create", &GeoModelBlueprintCreater::create);

    py::class_<GeoModelBlueprintCreater::Config>(bpc, "Config")
        .def(py::init<>())
        .def_readwrite("detectorSurfaces",
                       &GeoModelBlueprintCreater::Config::detectorSurfaces)
        .def_readwrite("kdtBinning",
                       &GeoModelBlueprintCreater::Config::kdtBinning);

    py::class_<GeoModelBlueprintCreater::Options>(bpc, "Options")
        .def(py::init<>())
        .def_readwrite("topEntry", &GeoModelBlueprintCreater::Options::topEntry)
        .def_readwrite("topBoundsOverride",
                       &GeoModelBlueprintCreater::Options::topBoundsOverride)
        .def_readwrite("table", &GeoModelBlueprintCreater::Options::table)
        .def_readwrite("dotGraph",
                       &GeoModelBlueprintCreater::Options::dotGraph);
  }

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
          return GeoModelDetectorElement::createDetectorElement<DiscSurface,
                                                                AnnulusBounds>(
              detElement->physicalVolume(), bounds, trafo,
              detElement->thickness());
        };

        return ITk::splitDiscModule(gctx, detElement, patterns, factory, name,
                                    *logger);
      },
      "gxtx"_a, "detElement"_a, "splitRanges"_a, "logLevel"_a = Logging::INFO);

  py::class_<GeoModelDetectorElementITk,
             std::shared_ptr<GeoModelDetectorElementITk>>(
      gm, "GeoModelDetectorElementITk")
      .def("surface", [](GeoModelDetectorElementITk& self) {
        return self.surface().getSharedPtr();
      });
  gm.def("convertToItk", &GeoModelDetectorElementITk::convertFromGeomodel);
}

}  // namespace ActsPython
