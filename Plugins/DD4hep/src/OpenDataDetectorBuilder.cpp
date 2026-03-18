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
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"
#include "ActsPlugins/Root/BlueprintBuilder.hpp"

#include <format>
#include <memory>
#include <optional>
#include <regex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DD4hep/DetType.h>

namespace {

const std::regex kPixelLayerFilter{"(?:PixelLayer|PixelEndcap[NP])(\\d)"};
const std::regex kShortStripLayerFilter{
    "(?:ShortStripLayer|ShortStripEndcap[NP])(\\d)"};
const std::regex kLongStripLayerFilter{
    "(?:LongStripLayer|LongStripEndcap[NP])(\\d)"};

static const Acts::ExtentEnvelope kBlueprintEnvelope =
    Acts::ExtentEnvelope::Zero()
        .set(Acts::AxisDirection::AxisZ, {20., 20.})
        .set(Acts::AxisDirection::AxisR, {0., 20.});

static const Acts::ExtentEnvelope kLayerEnvelope =
    Acts::ExtentEnvelope::Zero()
        .set(Acts::AxisDirection::AxisZ, {2., 2.})
        .set(Acts::AxisDirection::AxisR, {2., 2.});

template <typename Builder>
auto makeLayerCustomizer(Builder& builder, std::string det,
                         std::regex layerFilter) {
  using Element = typename Builder::Element;

  return [&builder, det = std::move(det), layerFilter = std::move(layerFilter)](
             const std::optional<Element>& elem,
             Acts::Experimental::LayerBlueprintNode& layer) {
    layer.setEnvelope(kLayerEnvelope);

    int layerIdx = 0;
    std::cmatch match;
    const std::string elemName =
        elem.has_value() ? std::string{builder.backend().nameOf(*elem)}
                         : layer.name();
    if (std::regex_search(elemName.c_str(), match, layerFilter) &&
        match.size() > 1) {
      layerIdx = std::stoi(match[1].str());
    } else {
      static const std::regex groupedLayerNameFilter{"layer(\\d+)"};
      if (std::regex_search(elemName.c_str(), match, groupedLayerNameFilter) &&
          match.size() > 1) {
        layerIdx = std::stoi(match[1].str());
      }
    }

    using SrfArrayNavPol = Acts::SurfaceArrayNavigationPolicy;
    using enum SrfArrayNavPol::LayerType;

    const auto& backend = builder.backend();
    SrfArrayNavPol::Config navCfg;
    const bool isBarrelLayer =
        elem.has_value()
            ? backend.isBarrel(*elem)
            : layer.layerType() ==
                  Acts::Experimental::LayerBlueprintNode::LayerType::Cylinder;
    if (isBarrelLayer) {
      navCfg.layerType = Cylinder;
      navCfg.bins = {backend.constant("{}_b{}_sf_b_phi", det, layerIdx),
                     backend.constant("{}_b_sf_b_z", det)};
    } else {
      navCfg.layerType = Disc;
      navCfg.bins = {backend.constant("{}_e_sf_b_r", det),
                     backend.constant("{}_e_sf_b_phi", det)};
    }

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

std::string pathFromElement(const ActsPlugins::TGeoBackend::Element& element) {
  std::vector<std::string> names;
  for (const auto* current = element.context.get(); current != nullptr;
       current = current->parent.get()) {
    if (current->node == nullptr) {
      continue;
    }
    if (current->node->GetName() != nullptr) {
      names.emplace_back(current->node->GetName());
    } else {
      names.emplace_back(current->node->GetVolume()->GetName());
    }
  }

  std::string path;
  for (auto it = names.rbegin(); it != names.rend(); ++it) {
    if (!path.empty()) {
      path += "|";
    }
    path += *it;
  }
  return path;
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
  cfg.sensitivePredicate =
      [nodeMap](const ActsPlugins::TGeoBackend::Element& element) {
        auto detElement = lookupDetElement(*nodeMap, element);
        return detElement.has_value() && detElement->volume().isSensitive();
      };
  cfg.barrelPredicate =
      [nodeMap](const ActsPlugins::TGeoBackend::Element& element) {
        auto detElement = lookupDetElement(*nodeMap, element);
        return detElement.has_value() &&
               dd4hep::DetType{detElement->typeFlag()}.is(
                   dd4hep::DetType::BARREL);
      };
  cfg.endcapPredicate =
      [nodeMap](const ActsPlugins::TGeoBackend::Element& element) {
        auto detElement = lookupDetElement(*nodeMap, element);
        return detElement.has_value() &&
               dd4hep::DetType{detElement->typeFlag()}.is(
                   dd4hep::DetType::ENDCAP);
      };
  cfg.trackerPredicate =
      [nodeMap](const ActsPlugins::TGeoBackend::Element& element) {
        auto detElement = lookupDetElement(*nodeMap, element);
        return detElement.has_value() &&
               dd4hep::DetType{detElement->typeFlag()}.is(
                   dd4hep::DetType::TRACKER);
      };
  cfg.beampipePredicate =
      [nodeMap](const ActsPlugins::TGeoBackend::Element& element) {
        auto detElement = lookupDetElement(*nodeMap, element);
        return detElement.has_value() &&
               dd4hep::DetType{detElement->typeFlag()}.is(
                   dd4hep::DetType::BEAMPIPE);
      };
  cfg.identifierProvider =
      [nodeMap](const ActsPlugins::TGeoBackend::Element& element) {
        auto detElement = lookupDetElement(*nodeMap, element);
        if (detElement.has_value()) {
          return static_cast<ActsPlugins::TGeoDetectorElement::Identifier>(
              detElement->volumeID());
        }
        return static_cast<ActsPlugins::TGeoDetectorElement::Identifier>(
            std::hash<std::string>{}(pathFromElement(element)));
      };
  cfg.constantProvider = [&detector](std::string_view name) {
    return detector.constant<int>(std::string{name});
  };
  return cfg;
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
          return std::format("layer{}", match[1].str());
        }
        current = builder.backend().parent(current);
      }
      return std::string{builder.backend().nameOf(elem)};
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
