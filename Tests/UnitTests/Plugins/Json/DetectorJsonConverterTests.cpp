// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/CylindricalContainerBuilder.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "ActsPlugins/Json/DetectorJsonConverter.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <fstream>
#include <memory>
#include <numbers>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

/// Helper method that allows to use the already existing testing
/// infrastructure with the new const-correct detector design
///
std::vector<std::shared_ptr<Surface>> unpackSurfaces(
    const std::vector<Surface*>& surfaces) {
  std::vector<std::shared_ptr<Surface>> uSurfaces;
  uSurfaces.reserve(surfaces.size());
  for (auto* s : surfaces) {
    uSurfaces.push_back(s->getSharedPtr());
  }
  return uSurfaces;
}

DetectorJsonConverter::Options detrayOptions() {
  // Detray format test - manipulate for detray
  DetectorVolumeJsonConverter::Options detrayOptions;
  detrayOptions.transformOptions.writeIdentity = true;
  detrayOptions.transformOptions.transpose = true;
  detrayOptions.surfaceOptions.transformOptions =
      detrayOptions.transformOptions;
  detrayOptions.portalOptions.surfaceOptions = detrayOptions.surfaceOptions;
  return DetectorJsonConverter::Options{detrayOptions};
}

GeometryContext tContext;
auto cGeometry = CylindricalTrackingGeometry(tContext);

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(SingleEmptyVolumeDetector) {
  auto portalGenerator = Experimental::defaultPortalGenerator();

  // Create a single cylindrical volume
  Transform3 nominal = Transform3::Identity();
  auto bounds = std::make_unique<CylinderVolumeBounds>(0., 50., 100.);

  auto volume = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume", nominal, std::move(bounds),
      Experimental::tryAllPortals());

  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes = {volume};

  Experimental::GeometryIdGenerator::Config generatorConfig;
  Experimental::GeometryIdGenerator generator(
      generatorConfig,
      getDefaultLogger("SequentialIdGenerator", Logging::VERBOSE));
  auto cache = generator.generateCache();
  for (auto& vol : volumes) {
    generator.assignGeometryId(cache, *vol);
  }

  auto detector = Experimental::Detector::makeShared(
      "Detector", volumes, Experimental::tryRootVolumes());

  auto jDetector = DetectorJsonConverter::toJson(tContext, *detector);

  std::ofstream out;
  out.open("single-empty-volume-detector.json");
  out << jDetector.dump(4);
  out.close();
}

BOOST_AUTO_TEST_CASE(SingleVolumeOneSurfaceDetector) {
  auto portalGenerator = Experimental::defaultPortalGenerator();

  // Create a single cylindrical volume
  Transform3 nominal = Transform3::Identity();
  auto bounds = std::make_unique<CylinderVolumeBounds>(0., 50., 100.);

  auto cylinderBounds = std::make_shared<CylinderBounds>(30., 90.);
  auto surface = Surface::makeShared<CylinderSurface>(nominal, cylinderBounds);

  auto volume = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume", nominal, std::move(bounds),
      {surface}, {}, Experimental::tryNoVolumes(),
      Experimental::tryAllPortalsAndSurfaces());

  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes = {volume};

  Experimental::GeometryIdGenerator::Config generatorConfig;
  Experimental::GeometryIdGenerator generator(
      generatorConfig,
      getDefaultLogger("SequentialIdGenerator", Logging::VERBOSE));
  auto cache = generator.generateCache();
  for (auto& vol : volumes) {
    generator.assignGeometryId(cache, *vol);
  }

  auto detector = Experimental::Detector::makeShared(
      "Detector", volumes, Experimental::tryRootVolumes());

  auto jDetector = DetectorJsonConverter::toJson(tContext, *detector);

  std::ofstream out;
  out.open("single-volume-one-surface-detector.json");
  out << jDetector.dump(4);
  out.close();

  out.open("single-volume-one-surface-detector-detray.json");
  out << DetectorJsonConverter::toJsonDetray(tContext, *detector,
                                             detrayOptions())
             .dump(4);
  out.close();
}

BOOST_AUTO_TEST_CASE(BeamPipeEndcapBarrelDetector) {
  // Detector store
  CylindricalTrackingGeometry::DetectorStore dStore;

  // Endcaps
  std::vector<std::shared_ptr<Experimental::IDetectorComponentBuilder>>
      endcapBuilders;
  for (auto [ie, ep] : enumerate(std::vector<double>({-710., 710.}))) {
    auto rSurfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 36., 0.125, 0.,
                                            55., ep, 2., 22u);

    auto endcapSurfaces =
        std::make_shared<Experimental::LayerStructureBuilder::SurfacesHolder>(
            unpackSurfaces(rSurfaces));
    // Configure the layer structure builder
    Experimental::LayerStructureBuilder::Config lsConfig;
    lsConfig.auxiliary = "*** Endcap with 22 surfaces ***";
    lsConfig.surfacesProvider = endcapSurfaces;
    lsConfig.binnings = {
        {DirectedProtoAxis(AxisDirection::AxisPhi, AxisBoundaryType::Closed,
                           -std::numbers::pi, std::numbers::pi, 22u),
         1u}};

    auto layerBuilder = std::make_shared<Experimental::LayerStructureBuilder>(
        lsConfig, getDefaultLogger("EndcapInteralsBuilder", Logging::VERBOSE));

    Experimental::VolumeStructureBuilder::Config shapeConfig;
    shapeConfig.boundValues = {18, 100, 10., std::numbers::pi, 0.};
    shapeConfig.transform =
        Transform3{Transform3::Identity()}.pretranslate(Vector3(0., 0., ep));
    shapeConfig.boundsType = VolumeBounds::BoundsType::eCylinder;

    auto shapeBuilder = std::make_shared<Experimental::VolumeStructureBuilder>(
        shapeConfig, getDefaultLogger("EndcapShapeBuilder", Logging::VERBOSE));

    Experimental::DetectorVolumeBuilder::Config dvCfg;
    dvCfg.name = "EndcapWithSurfaces_" + std::to_string(ie);
    dvCfg.externalsBuilder = shapeBuilder;
    dvCfg.internalsBuilder = layerBuilder;

    auto dvBuilder = std::make_shared<Experimental::DetectorVolumeBuilder>(
        dvCfg, getDefaultLogger("EndcapBuilder", Logging::VERBOSE));
    endcapBuilders.push_back(dvBuilder);
  }

  // Central barrel
  Experimental::VolumeStructureBuilder::Config innerShapeConfig;
  innerShapeConfig.boundValues = {18., 60., 700., std::numbers::pi, 0.};
  innerShapeConfig.boundsType = VolumeBounds::BoundsType::eCylinder;

  auto innerShapeBuilder =
      std::make_shared<Experimental::VolumeStructureBuilder>(
          innerShapeConfig,
          getDefaultLogger("InnerShapeBuilder", Logging::VERBOSE));

  Experimental::DetectorVolumeBuilder::Config ivCfg;
  ivCfg.name = "InnerBarrelGap";
  ivCfg.externalsBuilder = innerShapeBuilder;

  auto ivBuilder = std::make_shared<Experimental::DetectorVolumeBuilder>(
      ivCfg, getDefaultLogger("InnerBarrel", Logging::VERBOSE));

  ///  Barrel surfaces
  auto cSurfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.145, 72,
                                              3., 2., {32u, 14u});

  auto barrelSurfaces =
      std::make_shared<Experimental::LayerStructureBuilder::SurfacesHolder>(
          unpackSurfaces(cSurfaces));

  // Configure the layer structure builder
  Experimental::LayerStructureBuilder::Config lsConfig;
  lsConfig.auxiliary = "*** Barrel with 448 surfaces ***";
  lsConfig.surfacesProvider = barrelSurfaces;
  lsConfig.binnings = {
      {DirectedProtoAxis(AxisDirection::AxisZ, AxisBoundaryType::Bound, -480.,
                         480., 14u),
       1u},
      {DirectedProtoAxis(AxisDirection::AxisPhi, AxisBoundaryType::Closed,
                         -std::numbers::pi, std::numbers::pi, 32u),
       1u}};

  auto barrelBuilder = std::make_shared<Experimental::LayerStructureBuilder>(
      lsConfig, getDefaultLogger("BarrelInternalsBuilder", Logging::VERBOSE));

  Experimental::VolumeStructureBuilder::Config shapeConfig;
  shapeConfig.boundValues = {60., 80., 700., std::numbers::pi, 0.};
  shapeConfig.boundsType = VolumeBounds::BoundsType::eCylinder;

  auto shapeBuilder = std::make_shared<Experimental::VolumeStructureBuilder>(
      shapeConfig, getDefaultLogger("BarrelShapeBuilder", Logging::VERBOSE));

  Experimental::DetectorVolumeBuilder::Config dvCfg;
  dvCfg.name = "BarrelWithSurfaces";
  dvCfg.externalsBuilder = shapeBuilder;
  dvCfg.internalsBuilder = barrelBuilder;

  auto dvBuilder = std::make_shared<Experimental::DetectorVolumeBuilder>(
      dvCfg, getDefaultLogger("BarrelBuilder", Logging::VERBOSE));

  // Outer shape
  Experimental::VolumeStructureBuilder::Config outerShapeConfig;
  outerShapeConfig.boundValues = {80., 100., 700., std::numbers::pi, 0.};
  outerShapeConfig.boundsType = VolumeBounds::BoundsType::eCylinder;

  auto outerShapeBuilder =
      std::make_shared<Experimental::VolumeStructureBuilder>(
          outerShapeConfig,
          getDefaultLogger("OuterShapeBuilder", Logging::VERBOSE));

  Experimental::DetectorVolumeBuilder::Config ovCfg;
  ovCfg.name = "OuterBarrelGap";
  ovCfg.externalsBuilder = outerShapeBuilder;

  auto ovBuilder = std::make_shared<Experimental::DetectorVolumeBuilder>(
      ovCfg, getDefaultLogger("OuterBarrel", Logging::VERBOSE));

  // Build the combined barrel
  Experimental::CylindricalContainerBuilder::Config ccBarrelBuilderCfg;
  ccBarrelBuilderCfg.builders = {ivBuilder, dvBuilder, ovBuilder};
  ccBarrelBuilderCfg.binning = {AxisDirection::AxisR};

  auto ccBarrelBuilder =
      std::make_shared<Experimental::CylindricalContainerBuilder>(
          ccBarrelBuilderCfg,
          getDefaultLogger("BarrelBuilder", Logging::VERBOSE));

  // Builder the combined endcap barrel system
  Experimental::CylindricalContainerBuilder::Config ccBarrelEcBuilderCfg;
  ccBarrelEcBuilderCfg.builders = {endcapBuilders[0u], ccBarrelBuilder,
                                   endcapBuilders[1u]};
  ccBarrelEcBuilderCfg.binning = {AxisDirection::AxisZ};

  auto ccBarrelEndcapBuilder =
      std::make_shared<Experimental::CylindricalContainerBuilder>(
          ccBarrelEcBuilderCfg,
          getDefaultLogger("BarrelEndcapBuilder", Logging::VERBOSE));

  // Beam Pipe
  Experimental::VolumeStructureBuilder::Config bpShapeConfig;
  bpShapeConfig.boundValues = {0., 18., 720., std::numbers::pi, 0.};
  bpShapeConfig.boundsType = VolumeBounds::BoundsType::eCylinder;

  auto bpShapeBuilder = std::make_shared<Experimental::VolumeStructureBuilder>(
      bpShapeConfig,
      getDefaultLogger("BeamPipeShapeBuilder", Logging::VERBOSE));

  Experimental::DetectorVolumeBuilder::Config bpCfg;
  bpCfg.name = "BeamPipe";
  bpCfg.externalsBuilder = bpShapeBuilder;

  auto bpBuilder = std::make_shared<Experimental::DetectorVolumeBuilder>(
      bpCfg, getDefaultLogger("BeamPipe", Logging::VERBOSE));

  // Full detector
  Experimental::CylindricalContainerBuilder::Config detCompBuilderCfg;
  detCompBuilderCfg.builders = {bpBuilder, ccBarrelEndcapBuilder};
  detCompBuilderCfg.binning = {AxisDirection::AxisR};

  auto detCompBuilder =
      std::make_shared<Experimental::CylindricalContainerBuilder>(
          detCompBuilderCfg);

  auto gigConfig = Experimental::GeometryIdGenerator::Config();
  auto gig = std::make_shared<Experimental::GeometryIdGenerator>(gigConfig);

  Experimental::DetectorBuilder::Config detBuilderCfg;
  detBuilderCfg.name = "Detector";
  detBuilderCfg.builder = detCompBuilder;
  detBuilderCfg.geoIdGenerator = gig;

  auto detBuilder =
      std::make_shared<Experimental::DetectorBuilder>(detBuilderCfg);

  auto detector = detBuilder->construct(tContext);

  auto jDetector = DetectorJsonConverter::toJson(tContext, *detector);

  std::ofstream out;

  out.open("barrel-endcap-detector.json");
  out << jDetector.dump(4);
  out.close();

  auto in = std::ifstream("barrel-endcap-detector.json",
                          std::ifstream::in | std::ifstream::binary);

  BOOST_CHECK(in.good());
  nlohmann::json jDetectorIn;
  in >> jDetectorIn;
  in.close();

  auto detectorIn = DetectorJsonConverter::fromJson(tContext, jDetectorIn);

  BOOST_CHECK_EQUAL(detectorIn->name(), detector->name());

  auto jDetectorInOut = DetectorJsonConverter::toJson(tContext, *detectorIn);

  out.open("barrel-endcap-detector-closure.json");
  out << jDetectorInOut.dump(4);
  out.close();

  auto jDetectorDetray =
      DetectorJsonConverter::toJsonDetray(tContext, *detector, detrayOptions());

  out.open("barrel-endcap-detector-detray.json");
  out << jDetectorDetray.dump(4);
  out.close();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
