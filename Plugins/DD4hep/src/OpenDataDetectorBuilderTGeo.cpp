// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "ActsPlugins/DD4hep/OpenDataDetectorBuilder.hpp"
#include "ActsPlugins/Root/BlueprintBuilder.hpp"

#include <array>
#include <format>
#include <memory>
#include <optional>
#include <regex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>

#include <DD4hep/Detector.h>

#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"

namespace {

int oddTGeoConstant(std::string_view name) {
  static const std::unordered_map<std::string_view, int> kConstants = {
      {"pix_e_sf_b_r", 2},     {"pix_e_sf_b_phi", 36},  {"pix_b_sf_b_z", 14},
      {"pix_b0_sf_b_phi", 16}, {"pix_b1_sf_b_phi", 32}, {"pix_b2_sf_b_phi", 52},
      {"pix_b3_sf_b_phi", 78}, {"ss_e_sf_b_r", 3},      {"ss_e_sf_b_phi", 42},
      {"ss_b_sf_b_z", 21},     {"ss_b0_sf_b_phi", 40},  {"ss_b1_sf_b_phi", 56},
      {"ss_b2_sf_b_phi", 78},  {"ss_b3_sf_b_phi", 102}, {"ls_e_sf_b_r", 2},
      {"ls_e_sf_b_phi", 48},   {"ls_b_sf_b_z", 21},     {"ls_b0_sf_b_phi", 60},
      {"ls_b1_sf_b_phi", 80},
  };

  if (auto it = kConstants.find(name); it != kConstants.end()) {
    return it->second;
  }

  throw std::runtime_error(
      std::format("Unknown hard-coded TGeo ODD constant '{}'", name));
}

std::array<std::size_t, 2> oddTGeoBins(std::string_view det, bool isBarrelLayer,
                                       int layerIdx) {
  if (isBarrelLayer) {
    return {static_cast<std::size_t>(
                oddTGeoConstant(std::format("{}_b{}_sf_b_phi", det, layerIdx))),
            static_cast<std::size_t>(
                oddTGeoConstant(std::format("{}_b_sf_b_z", det)))};
  }

  return {static_cast<std::size_t>(
              oddTGeoConstant(std::format("{}_e_sf_b_r", det))),
          static_cast<std::size_t>(
              oddTGeoConstant(std::format("{}_e_sf_b_phi", det)))};
}

bool hasSensitiveMaterial(const ActsPlugins::TGeoBackend::Element& element) {
  if (element.context == nullptr || element.context->node == nullptr) {
    return false;
  }

  const auto* volume = element.context->node->GetVolume();
  if (volume == nullptr) {
    return false;
  }

  const auto* medium = volume->GetMedium();
  if (medium == nullptr || medium->GetMaterial() == nullptr ||
      medium->GetMaterial()->GetName() == nullptr) {
    return false;
  }

  return std::string_view{medium->GetMaterial()->GetName()}.find("Silicon") !=
         std::string_view::npos;
}

auto makeTGeoLayerCustomizer(ActsPlugins::BlueprintBuilder& builder,
                             std::string det, std::regex layerFilter) {
  return [&builder, det = std::move(det), layerFilter = std::move(layerFilter)](
             const std::optional<ActsPlugins::TGeoBackend::Element>& elem,
             Acts::Experimental::LayerBlueprintNode& layer) {
    layer.setEnvelope(ActsPlugins::DD4hep::detail::kLayerEnvelope);

    const std::string elemName =
        elem.has_value() ? std::string{builder.backend().nameOf(*elem)}
                         : layer.name();
    const int layerIdx =
        ActsPlugins::DD4hep::detail::layerIndexFromName(elemName, layerFilter);

    using SrfArrayNavPol = Acts::SurfaceArrayNavigationPolicy;
    using enum SrfArrayNavPol::LayerType;

    SrfArrayNavPol::Config navCfg;
    const bool isBarrelLayer =
        layer.layerType() ==
        Acts::Experimental::LayerBlueprintNode::LayerType::Cylinder;
    navCfg.layerType = isBarrelLayer ? Cylinder : Disc;

    const auto bins = oddTGeoBins(det, isBarrelLayer, layerIdx);
    navCfg.bins = {bins[0], bins[1]};

    layer.setNavigationPolicyFactory(Acts::NavigationPolicyFactory{}
                                         .add<Acts::CylinderNavigationPolicy>()
                                         .add<SrfArrayNavPol>(navCfg)
                                         .asUniquePtr());
  };
}

ActsPlugins::TGeoBackend::Config makeTGeoConfig(
    const dd4hep::Detector& detector) {
  ActsPlugins::TGeoBackend::Config cfg;
  cfg.root = detector.world().placement().ptr();
  cfg.lengthScale = Acts::UnitConstants::cm;
  cfg.sensitivePredicate = hasSensitiveMaterial;
  return cfg;
}

std::shared_ptr<Acts::Experimental::StaticBlueprintNode>
makeTemporaryTGeoBeampipeNode() {
  auto volumeBounds = std::make_shared<Acts::CylinderVolumeBounds>(
      0., 19. * Acts::UnitConstants::mm, 4. * Acts::UnitConstants::m);
  std::unique_ptr volume = std::make_unique<Acts::TrackingVolume>(
      Acts::Transform3::Identity(), volumeBounds, "BeamPipe");
  return std::make_shared<Acts::Experimental::StaticBlueprintNode>(
      std::move(volume));
}

void configureSubsystemNode(Acts::Experimental::ContainerBlueprintNode& node) {
  node.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
  node.setResizeStrategies(Acts::VolumeResizeStrategy::Gap,
                           Acts::VolumeResizeStrategy::Gap);
}

struct TGeoSubsystemSpec {
  std::string_view assembly;
  std::string_view det;
  std::string_view barrelName;
  std::string_view negativeEndcapName;
  std::string_view positiveEndcapName;
  const std::regex& layerFilter;
  const std::regex& barrelLayerFilter;
  const std::regex& negativeEndcapLayerFilter;
  const std::regex& positiveEndcapLayerFilter;
};

void addTGeoSubsystem(ActsPlugins::BlueprintBuilder& builder,
                      Acts::Experimental::BlueprintNode& parent,
                      const TGeoSubsystemSpec& spec) {
  const auto assemblyElement =
      builder.findDetElementByName(std::string{spec.assembly});
  if (!assemblyElement.has_value()) {
    throw std::runtime_error(
        std::format("Could not find assembly '{}'", spec.assembly));
  }
  const auto barrelElement = builder.findDetElementByName(
      *assemblyElement, std::string{spec.barrelName});
  if (!barrelElement.has_value()) {
    throw std::runtime_error(
        std::format("Could not find barrel '{}'", spec.barrelName));
  }
  const auto negativeEndcapElement = builder.findDetElementByName(
      *assemblyElement, std::string{spec.negativeEndcapName});
  if (!negativeEndcapElement.has_value()) {
    throw std::runtime_error(std::format("Could not find negative endcap '{}'",
                                         spec.negativeEndcapName));
  }
  const auto positiveEndcapElement = builder.findDetElementByName(
      *assemblyElement, std::string{spec.positiveEndcapName});
  if (!positiveEndcapElement.has_value()) {
    throw std::runtime_error(std::format("Could not find positive endcap '{}'",
                                         spec.positiveEndcapName));
  }

  auto subsystemNode =
      std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
          builder.backend().nameOf(*assemblyElement),
          Acts::AxisDirection::AxisZ);

  auto barrelNode = builder.layers()
                        .barrel()
                        .setSensorAxes("XYZ")
                        .setLayerFilter(spec.barrelLayerFilter)
                        .setContainer(*barrelElement)
                        .setContainerName(std::string{spec.barrelName})
                        .onLayer(makeTGeoLayerCustomizer(
                            builder, std::string{spec.det}, spec.layerFilter))
                        .build();
  configureSubsystemNode(*barrelNode);
  subsystemNode->addChild(std::move(barrelNode));

  auto negativeEndcapNode =
      builder.layers()
          .endcap()
          .setSensorAxes("XZY")
          .setLayerFilter(spec.negativeEndcapLayerFilter)
          .setContainer(*negativeEndcapElement)
          .setContainerName(std::string{spec.negativeEndcapName})
          .onLayer(makeTGeoLayerCustomizer(builder, std::string{spec.det},
                                           spec.layerFilter))
          .build();
  configureSubsystemNode(*negativeEndcapNode);
  subsystemNode->addChild(std::move(negativeEndcapNode));

  auto positiveEndcapNode =
      builder.layers()
          .endcap()
          .setSensorAxes("XZY")
          .setLayerFilter(spec.positiveEndcapLayerFilter)
          .setContainer(*positiveEndcapElement)
          .setContainerName(std::string{spec.positiveEndcapName})
          .onLayer(makeTGeoLayerCustomizer(builder, std::string{spec.det},
                                           spec.layerFilter))
          .build();
  configureSubsystemNode(*positiveEndcapNode);
  subsystemNode->addChild(std::move(positiveEndcapNode));

  parent.addChild(std::move(subsystemNode));
}

}  // namespace

namespace ActsPlugins::DD4hep {

std::unique_ptr<Acts::TrackingGeometry>
buildOpenDataDetectorBarrelEndcapViaTGeo(const dd4hep::Detector& detector,
                                         const Acts::GeometryContext& gctx,
                                         const Acts::Logger& logger) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using enum AxisDirection;

  ActsPlugins::BlueprintBuilder builder{makeTGeoConfig(detector),
                                        logger.cloneWithSuffix("TGeoBlpBld")};

  Blueprint::Config blueprintCfg;
  blueprintCfg.envelope = ActsPlugins::DD4hep::detail::kBlueprintEnvelope;
  Blueprint root{blueprintCfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  outer.addChild(makeTemporaryTGeoBeampipeNode());
  addTGeoSubsystem(
      builder, outer,
      {.assembly = "Pixels",
       .det = "pix",
       .barrelName = "PixelBarrel",
       .negativeEndcapName = "PixelEndcapN",
       .positiveEndcapName = "PixelEndcapP",
       .layerFilter = ActsPlugins::DD4hep::detail::kTGeoPixelLayerFilter,
       .barrelLayerFilter =
           ActsPlugins::DD4hep::detail::kTGeoPixelBarrelLayerFilter,
       .negativeEndcapLayerFilter =
           ActsPlugins::DD4hep::detail::kPixelNegativeEndcapLayerFilter,
       .positiveEndcapLayerFilter =
           ActsPlugins::DD4hep::detail::kPixelPositiveEndcapLayerFilter});
  addTGeoSubsystem(
      builder, outer,
      {.assembly = "ShortStrips",
       .det = "ss",
       .barrelName = "ShortStripBarrel",
       .negativeEndcapName = "ShortStripEndcapN",
       .positiveEndcapName = "ShortStripEndcapP",
       .layerFilter = ActsPlugins::DD4hep::detail::kTGeoShortStripLayerFilter,
       .barrelLayerFilter =
           ActsPlugins::DD4hep::detail::kTGeoShortStripBarrelLayerFilter,
       .negativeEndcapLayerFilter =
           ActsPlugins::DD4hep::detail::kShortStripNegativeEndcapLayerFilter,
       .positiveEndcapLayerFilter =
           ActsPlugins::DD4hep::detail::kShortStripPositiveEndcapLayerFilter});
  addTGeoSubsystem(
      builder, outer,
      {.assembly = "LongStrips",
       .det = "ls",
       .barrelName = "LongStripBarrel",
       .negativeEndcapName = "LongStripEndcapN",
       .positiveEndcapName = "LongStripEndcapP",
       .layerFilter = ActsPlugins::DD4hep::detail::kTGeoLongStripLayerFilter,
       .barrelLayerFilter =
           ActsPlugins::DD4hep::detail::kTGeoLongStripBarrelLayerFilter,
       .negativeEndcapLayerFilter =
           ActsPlugins::DD4hep::detail::kLongStripNegativeEndcapLayerFilter,
       .positiveEndcapLayerFilter =
           ActsPlugins::DD4hep::detail::kLongStripPositiveEndcapLayerFilter});

  return root.construct(BlueprintOptions{}, gctx, logger);
}

}  // namespace ActsPlugins::DD4hep
