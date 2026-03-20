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
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"

#include <format>
#include <memory>
#include <optional>
#include <regex>
#include <stdexcept>
#include <string>
#include <utility>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>

namespace ActsPlugins::DD4hep {

namespace {

auto makeLayerCustomizer(BlueprintBuilder& builder, std::string det,
                         std::regex layerFilter) {
  return [&builder, det = std::move(det), layerFilter = std::move(layerFilter)](
             const std::optional<dd4hep::DetElement>& elem,
             Acts::Experimental::LayerBlueprintNode& layer) {
    layer.setEnvelope(detail::kLayerEnvelope);

    const std::string elemName =
        elem.has_value() ? std::string{builder.backend().nameOf(*elem)}
                         : layer.name();
    const int layerIdx = detail::layerIndexFromName(elemName, layerFilter);

    using SrfArrayNavPol = Acts::SurfaceArrayNavigationPolicy;
    using enum SrfArrayNavPol::LayerType;

    SrfArrayNavPol::Config navCfg;
    const bool isBarrelLayer =
        layer.layerType() ==
        Acts::Experimental::LayerBlueprintNode::LayerType::Cylinder;
    if (isBarrelLayer) {
      navCfg.layerType = Cylinder;
      navCfg.bins = {
          builder.backend().constant("{}_b{}_sf_b_phi", det, layerIdx),
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

}  // namespace

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
  blueprintCfg.envelope = ActsPlugins::DD4hep::detail::kBlueprintEnvelope;
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

  addSubsystem("Pixels", "pix", ActsPlugins::DD4hep::detail::kPixelLayerFilter);
  addSubsystem("ShortStrips", "ss",
               ActsPlugins::DD4hep::detail::kShortStripLayerFilter);
  addSubsystem("LongStrips", "ls",
               ActsPlugins::DD4hep::detail::kLongStripLayerFilter);

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
  blueprintCfg.envelope = ActsPlugins::DD4hep::detail::kBlueprintEnvelope;
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

  addSubsystem("Pixels", "pix", ActsPlugins::DD4hep::detail::kPixelLayerFilter);
  addSubsystem("ShortStrips", "ss",
               ActsPlugins::DD4hep::detail::kShortStripLayerFilter);
  addSubsystem("LongStrips", "ls",
               ActsPlugins::DD4hep::detail::kLongStripLayerFilter);

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
  blueprintCfg.envelope = ActsPlugins::DD4hep::detail::kBlueprintEnvelope;
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

  addSubsystem("Pixels", "pix", ActsPlugins::DD4hep::detail::kPixelLayerFilter);
  addSubsystem("ShortStrips", "ss",
               ActsPlugins::DD4hep::detail::kShortStripLayerFilter);
  addSubsystem("LongStrips", "ls",
               ActsPlugins::DD4hep::detail::kLongStripLayerFilter);

  return root.construct(BlueprintOptions{}, gctx, logger);
}

}  // namespace ActsPlugins::DD4hep
