// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/OpenDataDetector.hpp"

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
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPlugins/Root/TGeoAxes.hpp"

#include <format>
#include <memory>
#include <optional>
#include <regex>
#include <string>
#include <utility>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>

namespace {

const std::regex kPixelLayerFilter{"(?:PixelLayer|PixelEndcap[NP])(\\d)"};
const std::regex kShortStripLayerFilter{
    "(?:ShortStripLayer|ShortStripEndcap[NP])(\\d)"};
const std::regex kLongStripLayerFilter{
    "(?:LongStripLayer|LongStripEndcap[NP])(\\d)"};

/// Returns a layer customizer callback that configures envelope, layer index
/// extraction, and navigation policy (phi/z bins for barrel, r/phi for endcap).
/// The callback is invoked once per layer during blueprint construction.
auto makeLayerCustomizer(ActsPlugins::DD4hep::BlueprintBuilder& builder,
                         const Acts::ExtentEnvelope& layerEnvelope,
                         auto constant, std::string det,
                         std::regex layerFilter) {
  return [&builder, layerEnvelope, constant, det = std::move(det),
          layerFilter = std::move(layerFilter)](
             const std::optional<dd4hep::DetElement>& elem,
             Acts::Experimental::LayerBlueprintNode& layer) {
    layer.setEnvelope(layerEnvelope);

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

    SrfArrayNavPol::Config navCfg;
    const bool isBarrelLayer =
        elem.has_value()
            ? builder.backend().isBarrel(*elem)
            : layer.layerType() ==
                  Acts::Experimental::LayerBlueprintNode::LayerType::Cylinder;
    if (isBarrelLayer) {
      navCfg.layerType = Cylinder;
      navCfg.bins = {constant("{}_b{}_sf_b_phi", det, layerIdx),
                     constant("{}_b_sf_b_z", det)};
    } else {
      navCfg.layerType = Disc;
      navCfg.bins = {constant("{}_e_sf_b_r", det),
                     constant("{}_e_sf_b_phi", det)};
    }

    layer.setNavigationPolicyFactory(Acts::NavigationPolicyFactory{}
                                         .add<Acts::CylinderNavigationPolicy>()
                                         .add<SrfArrayNavPol>(navCfg)
                                         .asUniquePtr());
  };
}

}  // namespace

namespace ActsExamples {
OpenDataDetector::OpenDataDetector(const Config& cfg,
                                   const Acts::GeometryContext& gctx)
    : DD4hepDetectorBase{cfg}, m_cfg{cfg} {
  ACTS_INFO("OpenDataDetector construct");
  construct(gctx);
}

auto OpenDataDetector::config() const -> const Config& {
  return m_cfg;
}

std::shared_ptr<ActsPlugins::DD4hepDetectorElement>
OpenDataDetector::defaultDetectorElementFactory(
    const dd4hep::DetElement& element, ActsPlugins::TGeoAxes axes,
    double scale) {
  return std::make_shared<ActsPlugins::DD4hepDetectorElement>(element, axes,
                                                              scale);
}

void OpenDataDetector::construct(const Acts::GeometryContext& gctx) {
  switch (m_cfg.constructionMethod) {
    case Config::ConstructionMethod::BarrelEndcap:
      constructBarrelEndcap(gctx);
      break;
    case Config::ConstructionMethod::DirectLayer:
      constructDirectLayer(gctx);
      break;
    case Config::ConstructionMethod::DirectLayerGrouped:
      constructDirectLayerGrouped(gctx);
      break;
  }
}

void OpenDataDetector::constructBarrelEndcap(
    const Acts::GeometryContext& gctx) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using namespace Acts::UnitLiterals;
  using enum AxisDirection;

  ActsPlugins::DD4hep::BlueprintBuilder builder{
      {
          .dd4hepDetector = &dd4hepDetector(),
          .lengthScale = Acts::UnitConstants::cm,
          .gctx = gctx,
      },
      logger().cloneWithSuffix("BlpBld")};

  Blueprint::Config blueprintCfg;
  blueprintCfg.envelope = m_cfg.blueprintEnvelope;
  Blueprint root{blueprintCfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  outer.addChild(builder.backend().makeBeampipe());

  auto constant = [this]<typename... Args>(std::format_string<Args...> fmt,
                                           Args&&... values) -> int {
    return dd4hepDetector().constant<int>(
        std::format(fmt, std::forward<Args>(values)...));
  };

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
        .onLayer(makeLayerCustomizer(builder, m_cfg.layerEnvelope, constant,
                                     std::move(det), layerFilter))
        .onContainer(
            [](const auto&,
               Acts::Experimental::CylinderContainerBlueprintNode& node) {
              node.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
              node.setResizeStrategies(VolumeResizeStrategy::Gap,
                                       VolumeResizeStrategy::Gap);
            })
        .addTo(outer);
  };

  addSubsystem("Pixels", "pix", kPixelLayerFilter);
  addSubsystem("ShortStrips", "ss", kShortStripLayerFilter);
  addSubsystem("LongStrips", "ls", kLongStripLayerFilter);

  m_trackingGeometry = root.construct(BlueprintOptions{}, gctx, logger());
}

void OpenDataDetector::constructDirectLayer(const Acts::GeometryContext& gctx) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using namespace Acts::UnitLiterals;
  using enum AxisDirection;

  ActsPlugins::DD4hep::BlueprintBuilder builder{
      {
          .dd4hepDetector = &dd4hepDetector(),
          .lengthScale = Acts::UnitConstants::cm,
          .gctx = gctx,
      },
      logger().cloneWithSuffix("BlpBld")};

  Blueprint::Config blueprintCfg;
  blueprintCfg.envelope = m_cfg.blueprintEnvelope;
  Blueprint root{blueprintCfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  outer.addChild(builder.backend().makeBeampipe());

  auto constant = [this]<typename... Args>(std::format_string<Args...> fmt,
                                           Args&&... values) -> int {
    return dd4hepDetector().constant<int>(
        std::format(fmt, std::forward<Args>(values)...));
  };

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

    auto layerCustomizer = makeLayerCustomizer(
        builder, m_cfg.layerEnvelope, constant, std::move(det), layerFilter);

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

  m_trackingGeometry = root.construct(BlueprintOptions{}, gctx, logger());
}

void OpenDataDetector::constructDirectLayerGrouped(
    const Acts::GeometryContext& gctx) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using namespace Acts::UnitLiterals;
  using enum AxisDirection;

  ActsPlugins::DD4hep::BlueprintBuilder builder{
      {
          .dd4hepDetector = &dd4hepDetector(),
          .lengthScale = Acts::UnitConstants::cm,
          .gctx = gctx,
      },
      logger().cloneWithSuffix("BlpBld")};

  Blueprint::Config blueprintCfg;
  blueprintCfg.envelope = m_cfg.blueprintEnvelope;
  Blueprint root{blueprintCfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  outer.addChild(builder.backend().makeBeampipe());

  auto constant = [this]<typename... Args>(std::format_string<Args...> fmt,
                                           Args&&... values) -> int {
    return dd4hepDetector().constant<int>(
        std::format(fmt, std::forward<Args>(values)...));
  };

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

    auto layerCustomizer = makeLayerCustomizer(
        builder, m_cfg.layerEnvelope, constant, std::move(det), layerFilter);

    // Walks the parent chain of a sensor element and returns the name of the
    // first ancestor whose name matches layerFilter, formatted as "layerN".
    // Falls back to the element's own name if no matching ancestor is found.
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

  m_trackingGeometry = root.construct(BlueprintOptions{}, gctx, logger());
}

}  // namespace ActsExamples
