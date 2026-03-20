// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/DD4hep/OpenDataDetectorBuilder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"
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
#include <utility>
#include <vector>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>

namespace {

const std::regex kPixelLayerFilter{"(?:PixelLayer|PixelEndcap[NP])(\\d)"};
const std::regex kShortStripLayerFilter{
    "(?:ShortStripLayer|ShortStripEndcap[NP])(\\d)"};
const std::regex kLongStripLayerFilter{
    "(?:LongStripLayer|LongStripEndcap[NP])(\\d)"};
const std::regex kPixelBarrelLayerFilter{"PixelLayer\\d"};
const std::regex kPixelNegativeEndcapLayerFilter{"PixelEndcapN\\d"};
const std::regex kPixelPositiveEndcapLayerFilter{"PixelEndcapP\\d"};
const std::regex kShortStripBarrelLayerFilter{"ShortStripLayer\\d"};
const std::regex kShortStripNegativeEndcapLayerFilter{"ShortStripEndcapN\\d"};
const std::regex kShortStripPositiveEndcapLayerFilter{"ShortStripEndcapP\\d"};
const std::regex kLongStripBarrelLayerFilter{"LongStripLayer\\d"};
const std::regex kLongStripNegativeEndcapLayerFilter{"LongStripEndcapN\\d"};
const std::regex kLongStripPositiveEndcapLayerFilter{"LongStripEndcapP\\d"};

static const Acts::ExtentEnvelope kBlueprintEnvelope =
    Acts::ExtentEnvelope::Zero()
        .set(Acts::AxisDirection::AxisZ, {20., 20.})
        .set(Acts::AxisDirection::AxisR, {0., 20.});

static const Acts::ExtentEnvelope kLayerEnvelope =
    Acts::ExtentEnvelope::Zero()
        .set(Acts::AxisDirection::AxisZ, {2., 2.})
        .set(Acts::AxisDirection::AxisR, {2., 2.});

int oddTGeoConstant(std::string_view name) {
  static const std::unordered_map<std::string_view, int> kConstants = {
      {"pix_e_sf_b_r", 2},    {"pix_e_sf_b_phi", 36},
      {"pix_b_sf_b_z", 14},   {"pix_b0_sf_b_phi", 16},
      {"pix_b1_sf_b_phi", 32},{"pix_b2_sf_b_phi", 52},
      {"pix_b3_sf_b_phi", 78},{"ss_e_sf_b_r", 3},
      {"ss_e_sf_b_phi", 42},  {"ss_b_sf_b_z", 21},
      {"ss_b0_sf_b_phi", 40}, {"ss_b1_sf_b_phi", 56},
      {"ss_b2_sf_b_phi", 78}, {"ss_b3_sf_b_phi", 102},
      {"ls_e_sf_b_r", 2},     {"ls_e_sf_b_phi", 48},
      {"ls_b_sf_b_z", 21},    {"ls_b0_sf_b_phi", 60},
      {"ls_b1_sf_b_phi", 80},
  };

  if (auto it = kConstants.find(name); it != kConstants.end()) {
    return it->second;
  }

  throw std::runtime_error(
      std::format("Unknown hard-coded TGeo ODD constant '{}'", name));
}

int layerIndexFromName(std::string_view elemName, const std::regex& layerFilter) {
  std::cmatch match;
  if (std::regex_search(elemName.begin(), elemName.end(), match, layerFilter) &&
      match.size() > 1) {
    return std::stoi(match[1].str());
  }

  static const std::regex groupedLayerNameFilter{"layer(\\d+)"};
  if (std::regex_search(elemName.begin(), elemName.end(), match,
                        groupedLayerNameFilter) &&
      match.size() > 1) {
    return std::stoi(match[1].str());
  }

  return 0;
}

std::array<std::size_t, 2> oddTGeoBins(std::string_view det, bool isBarrelLayer,
                                       int layerIdx) {
  if (isBarrelLayer) {
    return {static_cast<std::size_t>(
                oddTGeoConstant(std::format("{}_b{}_sf_b_phi", det, layerIdx))),
            static_cast<std::size_t>(oddTGeoConstant(
                std::format("{}_b_sf_b_z", det)))};
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

auto makeLayerCustomizer(ActsPlugins::DD4hep::BlueprintBuilder& builder,
                         std::string det, std::regex layerFilter) {
  return [&builder, det = std::move(det), layerFilter = std::move(layerFilter)](
             const std::optional<dd4hep::DetElement>& elem,
             Acts::Experimental::LayerBlueprintNode& layer) {
    layer.setEnvelope(kLayerEnvelope);

    const std::string elemName =
        elem.has_value() ? std::string{builder.backend().nameOf(*elem)}
                         : layer.name();
    const int layerIdx = layerIndexFromName(elemName, layerFilter);

    using SrfArrayNavPol = Acts::SurfaceArrayNavigationPolicy;
    using enum SrfArrayNavPol::LayerType;

    SrfArrayNavPol::Config navCfg;
    const bool isBarrelLayer =
        layer.layerType() ==
        Acts::Experimental::LayerBlueprintNode::LayerType::Cylinder;
    if (isBarrelLayer) {
      navCfg.layerType = Cylinder;
      navCfg.bins = {builder.backend().constant("{}_b{}_sf_b_phi", det, layerIdx),
                     builder.backend().constant("{}_b_sf_b_z", det)};
    } else {
      navCfg.layerType = Disc;
      navCfg.bins = {builder.backend().constant("{}_e_sf_b_r", det),
                     builder.backend().constant("{}_e_sf_b_phi", det)};
    }

    layer.setNavigationPolicyFactory(Acts::NavigationPolicyFactory{}
                                         .add<Acts::CylinderNavigationPolicy>()
                                         .add<SrfArrayNavPol>(navCfg)
                                         .asUniquePtr());
  };
}

auto makeTGeoLayerCustomizer(ActsPlugins::BlueprintBuilder& builder,
                             std::string det, std::regex layerFilter) {
  return [&builder, det = std::move(det), layerFilter = std::move(layerFilter)](
             const std::optional<ActsPlugins::TGeoBackend::Element>& elem,
             Acts::Experimental::LayerBlueprintNode& layer) {
    layer.setEnvelope(kLayerEnvelope);

    const std::string elemName =
        elem.has_value() ? std::string{builder.backend().nameOf(*elem)}
                         : layer.name();
    const int layerIdx = layerIndexFromName(elemName, layerFilter);

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

using NodeToDetElementMap =
    std::unordered_map<const TGeoNode*, dd4hep::DetElement>;

void collectNodeToDetElementMap(const dd4hep::DetElement& detElement,
                                NodeToDetElementMap& map) {
  if (detElement.isValid() && detElement.placement().isValid()) {
    map.try_emplace(detElement.placement().ptr(), detElement);
  }

  for (const auto& [name, child] : detElement.children()) {
    (void)name;
    collectNodeToDetElementMap(child, map);
  }
}

std::shared_ptr<const NodeToDetElementMap> buildNodeToDetElementMap(
    const dd4hep::Detector& detector) {
  auto map = std::make_shared<NodeToDetElementMap>();
  collectNodeToDetElementMap(detector.world(), *map);
  return map;
}

std::optional<dd4hep::DetElement> lookupDetElement(
    const NodeToDetElementMap& map, const ActsPlugins::TGeoBackend::Element& e) {
  if (e.context == nullptr || e.context->node == nullptr) {
    return std::nullopt;
  }
  if (auto it = map.find(e.context->node); it != map.end()) {
    return it->second;
  }
  return std::nullopt;
}

ActsPlugins::TGeoBackend::Config makeTGeoConfigFromDD4hep(
    const dd4hep::Detector& detector) {
  auto nodeMap = buildNodeToDetElementMap(detector);

  ActsPlugins::TGeoBackend::Config cfg;
  cfg.root = detector.world().placement().ptr();
  cfg.lengthScale = Acts::UnitConstants::cm;
  cfg.nameProvider = [nodeMap](const ActsPlugins::TGeoBackend::Element& element)
      -> std::string {
    if (auto detElement = lookupDetElement(*nodeMap, element);
        detElement.has_value()) {
      return detElement->name();
    }

    if (element.context != nullptr && element.context->node != nullptr &&
        element.context->node->GetName() != nullptr) {
      return element.context->node->GetName();
    }
    return {};
  };
  cfg.sensitivePredicate = hasSensitiveMaterial;
  return cfg;
}

std::shared_ptr<Acts::Experimental::StaticBlueprintNode>
makeTemporaryTGeoBeampipeNode(const dd4hep::Detector& detector) {
  (void)detector;
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

void addTGeoSubsystem(
    ActsPlugins::BlueprintBuilder& builder,
    Acts::Experimental::BlueprintNode& parent,
    const TGeoSubsystemSpec& spec) {
  const auto assemblyElement = builder.findDetElementByName(std::string{spec.assembly});
  if (!assemblyElement.has_value()) {
    throw std::runtime_error(
        std::format("Could not find assembly '{}'", spec.assembly));
  }
  const auto barrelElement =
      builder.findDetElementByName(std::string{spec.barrelName});
  if (!barrelElement.has_value()) {
    throw std::runtime_error(
        std::format("Could not find barrel '{}'", spec.barrelName));
  }
  const auto negativeEndcapElement =
      builder.findDetElementByName(std::string{spec.negativeEndcapName});
  if (!negativeEndcapElement.has_value()) {
    throw std::runtime_error(std::format("Could not find negative endcap '{}'",
                                         spec.negativeEndcapName));
  }
  const auto positiveEndcapElement =
      builder.findDetElementByName(std::string{spec.positiveEndcapName});
  if (!positiveEndcapElement.has_value()) {
    throw std::runtime_error(std::format("Could not find positive endcap '{}'",
                                         spec.positiveEndcapName));
  }

  auto subsystemNode =
      std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
          builder.backend().nameOf(*assemblyElement), Acts::AxisDirection::AxisZ);

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

  auto negativeEndcapNode = builder.layers()
                                .endcap()
                                .setSensorAxes("XZY")
                                .setLayerFilter(spec.negativeEndcapLayerFilter)
                                .setContainer(*negativeEndcapElement)
                                .setContainerName(
                                    std::string{spec.negativeEndcapName})
                                .onLayer(makeTGeoLayerCustomizer(
                                    builder, std::string{spec.det},
                                    spec.layerFilter))
                                .build();
  configureSubsystemNode(*negativeEndcapNode);
  subsystemNode->addChild(std::move(negativeEndcapNode));

  auto positiveEndcapNode = builder.layers()
                                .endcap()
                                .setSensorAxes("XZY")
                                .setLayerFilter(spec.positiveEndcapLayerFilter)
                                .setContainer(*positiveEndcapElement)
                                .setContainerName(
                                    std::string{spec.positiveEndcapName})
                                .onLayer(makeTGeoLayerCustomizer(
                                    builder, std::string{spec.det},
                                    spec.layerFilter))
                                .build();
  configureSubsystemNode(*positiveEndcapNode);
  subsystemNode->addChild(std::move(positiveEndcapNode));

  parent.addChild(std::move(subsystemNode));
}

}  // namespace

namespace ActsPlugins::DD4hep {

std::unique_ptr<Acts::TrackingGeometry> buildOpenDataDetectorBarrelEndcap(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
    const Acts::Logger& logger) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using enum AxisDirection;

  BlueprintBuilder builder{{
                               .dd4hepDetector = &detector,
                               .lengthScale = Acts::UnitConstants::cm,
                               .gctx = gctx,
                           },
                           logger.cloneWithSuffix("BlpBld")};

  Blueprint::Config blueprintCfg;
  blueprintCfg.envelope = kBlueprintEnvelope;
  Blueprint root{blueprintCfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  outer.addChild(builder.backend().makeBeampipe());

  auto addSubsystem = [&](std::string assembly, std::string det,
                          const std::regex& layerFilter) {
    const auto assemblyElement = builder.findDetElementByName(assembly);
    if (!assemblyElement.has_value()) {
      throw std::runtime_error(
          std::format("Could not find assembly '{}'", assembly));
    }

    builder.barrelEndcap()
        .setAssembly(*assemblyElement)
        .setSensorAxes("XYZ", "XZY")
        .setLayerFilter(layerFilter)
        .onLayer(makeLayerCustomizer(builder, std::move(det), layerFilter))
        .onContainer(
            [](const auto&, Acts::Experimental::ContainerBlueprintNode& node) {
              node.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
              node.setResizeStrategies(VolumeResizeStrategy::Gap,
                                       VolumeResizeStrategy::Gap);
            })
        .addTo(outer);
  };

  addSubsystem("Pixels", "pix", kPixelLayerFilter);
  addSubsystem("ShortStrips", "ss", kShortStripLayerFilter);
  addSubsystem("LongStrips", "ls", kLongStripLayerFilter);

  return root.construct(BlueprintOptions{}, gctx, logger);
}

std::unique_ptr<Acts::TrackingGeometry>
buildOpenDataDetectorBarrelEndcapViaTGeo(const dd4hep::Detector& detector,
                                         const Acts::GeometryContext& gctx,
                                         const Acts::Logger& logger) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using enum AxisDirection;

  ActsPlugins::BlueprintBuilder builder{makeTGeoConfigFromDD4hep(detector),
                                        logger.cloneWithSuffix("TGeoBlpBld")};

  Blueprint::Config blueprintCfg;
  blueprintCfg.envelope = kBlueprintEnvelope;
  Blueprint root{blueprintCfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  outer.addChild(makeTemporaryTGeoBeampipeNode(detector));
  addTGeoSubsystem(builder, outer,
                   {.assembly = "Pixels",
                    .det = "pix",
                    .barrelName = "PixelBarrel",
                    .negativeEndcapName = "PixelEndcapN",
                    .positiveEndcapName = "PixelEndcapP",
                    .layerFilter = kPixelLayerFilter,
                    .barrelLayerFilter = kPixelBarrelLayerFilter,
                    .negativeEndcapLayerFilter = kPixelNegativeEndcapLayerFilter,
                    .positiveEndcapLayerFilter = kPixelPositiveEndcapLayerFilter});
  addTGeoSubsystem(
      builder, outer,
      {.assembly = "ShortStrips",
       .det = "ss",
       .barrelName = "ShortStripBarrel",
       .negativeEndcapName = "ShortStripEndcapN",
       .positiveEndcapName = "ShortStripEndcapP",
       .layerFilter = kShortStripLayerFilter,
       .barrelLayerFilter = kShortStripBarrelLayerFilter,
       .negativeEndcapLayerFilter = kShortStripNegativeEndcapLayerFilter,
       .positiveEndcapLayerFilter = kShortStripPositiveEndcapLayerFilter});
  addTGeoSubsystem(
      builder, outer,
      {.assembly = "LongStrips",
       .det = "ls",
       .barrelName = "LongStripBarrel",
       .negativeEndcapName = "LongStripEndcapN",
       .positiveEndcapName = "LongStripEndcapP",
       .layerFilter = kLongStripLayerFilter,
       .barrelLayerFilter = kLongStripBarrelLayerFilter,
       .negativeEndcapLayerFilter = kLongStripNegativeEndcapLayerFilter,
       .positiveEndcapLayerFilter = kLongStripPositiveEndcapLayerFilter});

  return root.construct(BlueprintOptions{}, gctx, logger);
}

std::unique_ptr<Acts::TrackingGeometry> buildOpenDataDetectorDirectLayer(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
    const Acts::Logger& logger) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using enum AxisDirection;

  BlueprintBuilder builder{{
                               .dd4hepDetector = &detector,
                               .lengthScale = Acts::UnitConstants::cm,
                               .gctx = gctx,
                           },
                           logger.cloneWithSuffix("BlpBld")};

  Blueprint::Config blueprintCfg;
  blueprintCfg.envelope = kBlueprintEnvelope;
  Blueprint root{blueprintCfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  outer.addChild(builder.backend().makeBeampipe());

  auto addSubsystem = [&](std::string assembly, std::string det,
                          const std::regex& layerFilter) {
    const auto assemblyElement = builder.findDetElementByName(assembly);
    if (!assemblyElement.has_value()) {
      throw std::runtime_error(
          std::format("Could not find assembly '{}'", assembly));
    }

    auto barrels = builder.findBarrelElements(*assemblyElement);
    auto endcaps = builder.findEndcapElements(*assemblyElement);

    const std::string assemblyName{builder.backend().nameOf(*assemblyElement)};
    auto containerNode =
        std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
            assemblyName, Acts::AxisDirection::AxisZ);

    auto layerCustomizer =
        makeLayerCustomizer(builder, std::move(det), layerFilter);

    for (const auto& barrel : barrels) {
      auto barrelNode = builder.layers()
                            .barrel()
                            .setSensorAxes("XYZ")
                            .setLayerFilter(layerFilter)
                            .setContainer(barrel)
                            .onLayer(layerCustomizer)
                            .build();
      barrelNode->setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
      barrelNode->setResizeStrategies(VolumeResizeStrategy::Gap,
                                      VolumeResizeStrategy::Gap);
      containerNode->addChild(std::move(barrelNode));
    }

    for (const auto& endcap : endcaps) {
      auto endcapNode = builder.layers()
                            .endcap()
                            .setSensorAxes("XZY")
                            .setLayerFilter(layerFilter)
                            .setContainer(endcap)
                            .onLayer(layerCustomizer)
                            .build();
      endcapNode->setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
      endcapNode->setResizeStrategies(VolumeResizeStrategy::Gap,
                                      VolumeResizeStrategy::Gap);
      containerNode->addChild(std::move(endcapNode));
    }

    outer.addChild(std::move(containerNode));
  };

  addSubsystem("Pixels", "pix", kPixelLayerFilter);
  addSubsystem("ShortStrips", "ss", kShortStripLayerFilter);
  addSubsystem("LongStrips", "ls", kLongStripLayerFilter);

  return root.construct(BlueprintOptions{}, gctx, logger);
}

std::unique_ptr<Acts::TrackingGeometry> buildOpenDataDetectorDirectLayerGrouped(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
    const Acts::Logger& logger) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using enum AxisDirection;

  BlueprintBuilder builder{{
                               .dd4hepDetector = &detector,
                               .lengthScale = Acts::UnitConstants::cm,
                               .gctx = gctx,
                           },
                           logger.cloneWithSuffix("BlpBld")};

  Blueprint::Config blueprintCfg;
  blueprintCfg.envelope = kBlueprintEnvelope;
  Blueprint root{blueprintCfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  outer.addChild(builder.backend().makeBeampipe());

  auto addSubsystem = [&](std::string assembly, std::string det,
                          const std::regex& layerFilter) {
    const auto assemblyElement = builder.findDetElementByName(assembly);
    if (!assemblyElement.has_value()) {
      throw std::runtime_error(
          std::format("Could not find assembly '{}'", assembly));
    }

    auto barrels = builder.findBarrelElements(*assemblyElement);
    auto endcaps = builder.findEndcapElements(*assemblyElement);

    const std::string assemblyName{builder.backend().nameOf(*assemblyElement)};
    auto containerNode =
        std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
            assemblyName, Acts::AxisDirection::AxisZ);

    auto layerCustomizer =
        makeLayerCustomizer(builder, std::move(det), layerFilter);

    auto sensorToLayerKey = [&](const dd4hep::DetElement& elem) -> std::string {
      auto current = elem;
      const auto world = builder.backend().world();
      while (!(current == world)) {
        std::cmatch match;
        const std::string name{builder.backend().nameOf(current)};
        if (std::regex_search(name.c_str(), match, layerFilter) &&
            match.size() > 1) {
          return builder.getPathToElementName(current);
        }
        current = builder.backend().parent(current);
      }
      return builder.getPathToElementName(elem);
    };

    for (const auto& barrel : barrels) {
      auto sensors = builder.resolveSensitives(barrel);
      auto barrelNode = builder.layersFromSensors()
                            .barrel()
                            .setSensorAxes("XYZ")
                            .setSensors(std::move(sensors))
                            .setContainerName(builder.backend().nameOf(barrel))
                            .groupBy(sensorToLayerKey)
                            .onLayer(layerCustomizer)
                            .build();
      barrelNode->setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
      barrelNode->setResizeStrategies(VolumeResizeStrategy::Gap,
                                      VolumeResizeStrategy::Gap);
      containerNode->addChild(std::move(barrelNode));
    }

    for (const auto& endcap : endcaps) {
      auto sensors = builder.resolveSensitives(endcap);
      auto endcapNode = builder.layersFromSensors()
                            .endcap()
                            .setSensorAxes("XZY")
                            .setSensors(std::move(sensors))
                            .setContainerName(builder.backend().nameOf(endcap))
                            .groupBy(sensorToLayerKey)
                            .onLayer(layerCustomizer)
                            .build();
      endcapNode->setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
      endcapNode->setResizeStrategies(VolumeResizeStrategy::Gap,
                                      VolumeResizeStrategy::Gap);
      containerNode->addChild(std::move(endcapNode));
    }

    outer.addChild(std::move(containerNode));
  };

  addSubsystem("Pixels", "pix", kPixelLayerFilter);
  addSubsystem("ShortStrips", "ss", kShortStripLayerFilter);
  addSubsystem("LongStrips", "ls", kLongStripLayerFilter);

  return root.construct(BlueprintOptions{}, gctx, logger);
}

}  // namespace ActsPlugins::DD4hep
