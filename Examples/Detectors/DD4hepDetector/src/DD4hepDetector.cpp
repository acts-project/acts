// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "Acts/Detector/detail/BlueprintHelper.hpp"
#include "Acts/Detector/CylindricalContainerBuilder.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/DD4hep/DD4hepBlueprint.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorSurfaceFactory.hpp"
#include "Acts/Plugins/DD4hep/DD4hepFieldAdapter.hpp"
#include "Acts/Plugins/DD4hep/DD4hepLayerStructure.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Detector.h>
#include <DD4hep/Fields.h>
#include <boost/program_options.hpp>

namespace ActsExamples {
namespace DD4hep {

DD4hepDetector::DD4hepDetector() = default;

DD4hepDetector::DD4hepDetector(
    std::shared_ptr<DD4hepGeometryService> _geometryService)
    : geometryService(std::move(_geometryService)) {}

DD4hepDetector::~DD4hepDetector() = default;

auto DD4hepDetector::finalize(
    ActsExamples::DD4hep::DD4hepGeometryService::Config config,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  Acts::GeometryContext dd4HepContext;
  config.matDecorator = std::move(mdecorator);
  geometryService =
      std::make_shared<ActsExamples::DD4hep::DD4hepGeometryService>(config);
  TrackingGeometryPtr dd4tGeometry =
      geometryService->trackingGeometry(dd4HepContext);
  if (!dd4tGeometry) {
    throw std::runtime_error{
        "Did not receive tracking geometry from DD4hep geometry service"};
  }
  ContextDecorators dd4ContextDecorators = {};
  // return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(dd4tGeometry), std::move(dd4ContextDecorators));
}

auto DD4hepDetector::finalize(const Acts::GeometryContext& gctx,
              const std::vector<std::string>& compactFiles,
              std::shared_ptr<const Acts::IMaterialDecorator> /*unused*/)
    -> std::tuple<DetectorPtr, DetectorStore, ContextDecorators> {

  // Create the geometry seervice and the world from DD4hep
  auto lcdd = &(dd4hep::Detector::getInstance());
  for (const auto& compactFile : compactFiles) {
    lcdd->fromCompact(compactFile);
  }
  lcdd->volumeManager();
  lcdd->apply("DD4hepVolumeManager", 0, nullptr);
  auto world = lcdd->world();

  // Surface factory
  auto surfaceFactory = std::make_shared<Acts::DD4hepDetectorSurfaceFactory>(
      Acts::getDefaultLogger("DD4hepDetectorSurfaceFactory",
                             Acts::Logging::VERBOSE));

  // Layer structure geneterator
  auto layerStructure =
      std::make_shared<Acts::Experimental::DD4hepLayerStructure>(
          std::move(surfaceFactory),
          Acts::getDefaultLogger("DD4hepLayerStructure",
                                 Acts::Logging::VERBOSE));

  Acts::Experimental::DD4hepBlueprint::Config bpCfg{layerStructure};
  Acts::Experimental::DD4hepBlueprint::Cache bpCache;

  Acts::Experimental::DD4hepBlueprint blueprintCreator(
      bpCfg, Acts::getDefaultLogger("DD4hepBlueprint", Acts::Logging::VERBOSE));
  auto dd4hepBlueprint = blueprintCreator.create(bpCache, gctx, world);

  // Complete the gaps
  Acts::Experimental::detail::BlueprintHelper::fillGaps(*dd4hepBlueprint);

  // Create a Cylindrical detector builder from this blueprint
  auto dd4hepDtectorBuilder =
      std::make_shared<Acts::Experimental::CylindricalContainerBuilder>(
          *dd4hepBlueprint, Acts::Logging::VERBOSE);

  // Detector builder
  Acts::Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxiliary =
      "*** Test : auto generated cylindrical detector builder  ***";
  dCfg.name = "Cylindrical detector from blueprint";
  dCfg.builder = dd4hepDtectorBuilder;
  dCfg.geoIdGenerator = dd4hepBlueprint->geoIdGenerator;

  auto dd4hepDetector = Acts::Experimental::DetectorBuilder(dCfg).construct(gctx);

  // Return what you have
  ContextDecorators dd4ContextDecorators = {};
  return std::tie(dd4hepDetector, bpCache.dd4hepStore, dd4ContextDecorators);
}

std::shared_ptr<Acts::DD4hepFieldAdapter> DD4hepDetector::field() const {
  const auto& detector = geometryService->detector();

  return std::make_shared<Acts::DD4hepFieldAdapter>(detector.field());
}

}  // namespace DD4hep
}  // namespace ActsExamples
