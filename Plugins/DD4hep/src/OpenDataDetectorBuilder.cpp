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
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
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

// Placeholder bin counts for the (deferred) proto material grids attached to
// each layer / the beampipe. The axis *range* is resolved later from the
// actual surface bounds during Blueprint::construct; only the bin count is
// fixed here. TODO: derive these from the DD4hep XML (e.g. the same
// layer_material_*_bin* constants the Gen1 conversion reads) instead of a
// fixed placeholder.
constexpr std::size_t kMatPhiBins = 10;
constexpr std::size_t kMatZBins = 10;
constexpr std::size_t kMatRBins = 10;

auto makeLayerCustomizer(const BlueprintBuilder& builder, std::string det,
                         std::regex layerFilter,
                         Acts::CylinderVolumeBounds::Face barrelMaterialFace) {
  return [&builder, det = std::move(det), layerFilter = std::move(layerFilter),
          barrelMaterialFace](const std::optional<dd4hep::DetElement>& elem,
                              Acts::detail::LayerNodePtr layer)
             -> Acts::detail::BlueprintNodePtr {
    layer->setEnvelope(detail::kLayerEnvelope);

    const std::string elemName =
        elem.has_value() ? std::string{builder.backend().nameOf(*elem)}
                         : layer->name();
    const int layerIdx = detail::layerIndexFromName(elemName, layerFilter);

    using SrfArrayNavPol = Acts::SurfaceArrayNavigationPolicy;
    using enum SrfArrayNavPol::LayerType;

    SrfArrayNavPol::Config navCfg;
    navCfg.envelope = detail::kLayerEnvelope;

    using enum Acts::CylinderVolumeBounds::Face;
    auto matNode = std::make_shared<Acts::MaterialDesignatorBlueprintNode>(
        layer->name() + "_mat");

    if (layer->layerType() == Acts::LayerBlueprintNode::LayerType::Cylinder) {
      // Barrel layer
      navCfg.layerType = Cylinder;
      navCfg.bins = {
          builder.backend().constant("{}_b{}_sf_b_phi", det, layerIdx),
          builder.backend().constant("{}_b_sf_b_z", det)};

      // Barrel: thin cylindrical shell -> one mantle face carries material.
      // Which face (inner vs. outer) is a per-subsystem convention (Pixel
      // uses outer, ShortStrips/LongStrips use inner) chosen so that no two
      // radially-stacked layers independently claim the same fused portal.
      matNode->configureFace(barrelMaterialFace,
                             Acts::AxisSpec::DeferredEquidistant(
                                 kMatPhiBins, Acts::AxisDirection::AxisRPhi),
                             Acts::AxisSpec::DeferredEquidistant(
                                 kMatZBins, Acts::AxisDirection::AxisZ));
    } else {
      // Endcap layer
      navCfg.layerType = Disc;
      navCfg.bins = {builder.backend().constant("{}_e_sf_b_r", det),
                     builder.backend().constant("{}_e_sf_b_phi", det)};

      // Endcap: thin disc "pancake" -> both flat faces carry material
      matNode->configureFace(NegativeDisc,
                             Acts::AxisSpec::DeferredEquidistant(
                                 kMatRBins, Acts::AxisDirection::AxisR),
                             Acts::AxisSpec::DeferredEquidistant(
                                 kMatPhiBins, Acts::AxisDirection::AxisPhi));
      matNode->configureFace(PositiveDisc,
                             Acts::AxisSpec::DeferredEquidistant(
                                 kMatRBins, Acts::AxisDirection::AxisR),
                             Acts::AxisSpec::DeferredEquidistant(
                                 kMatPhiBins, Acts::AxisDirection::AxisPhi));
    }

    layer->setNavigationPolicyFactory(Acts::NavigationPolicyFactory{}
                                          .add<Acts::CylinderNavigationPolicy>()
                                          .add<SrfArrayNavPol>(navCfg)
                                          .asUniquePtr());

    matNode->addChild(std::move(layer));
    return matNode;
  };
}

void addDirectLayerSubsystem(
    const BlueprintBuilder& builder, Acts::ContainerBlueprintNode& outer,
    std::string assembly, std::string det, const std::regex& layerFilter,
    Acts::CylinderVolumeBounds::Face barrelMaterialFace) {
  const auto assemblyElement = builder.findDetElementByName(assembly);
  if (!assemblyElement.has_value()) {
    throw std::runtime_error(
        std::format("Could not find assembly '{}'", assembly));
  }

  auto barrels = builder.findBarrelElements(*assemblyElement);
  auto endcaps = builder.findEndcapElements(*assemblyElement);

  const std::string assemblyName{builder.backend().nameOf(*assemblyElement)};
  auto containerNode = std::make_shared<Acts::CylinderContainerBlueprintNode>(
      assemblyName, Acts::AxisDirection::AxisZ);

  auto layerCustomizer = makeLayerCustomizer(builder, std::move(det),
                                             layerFilter, barrelMaterialFace);

  auto addLayerChildren = [&](const auto& elements, auto makeNode) {
    for (const auto& element : elements) {
      auto node = makeNode(element);
      node->setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
      node->setResizeStrategies(Acts::VolumeResizeStrategy::Gap,
                                Acts::VolumeResizeStrategy::Gap);
      containerNode->addChild(std::move(node));
    }
  };

  addLayerChildren(barrels, [&](const auto& barrel) {
    return builder.layers()
        .barrel()
        .setSensorAxes("XYZ")
        .setLayerFilter(layerFilter)
        .setContainer(barrel)
        .onLayer(layerCustomizer)
        .build();
  });

  addLayerChildren(endcaps, [&](const auto& endcap) {
    return builder.layers()
        .endcap()
        .setSensorAxes("XZY")
        .setLayerFilter(layerFilter)
        .setContainer(endcap)
        .onLayer(layerCustomizer)
        .build();
  });

  outer.addChild(std::move(containerNode));
}

void addBarrelEndcapSubsystem(
    const BlueprintBuilder& builder, Acts::ContainerBlueprintNode& outer,
    std::string assembly, std::string det, const std::regex& layerFilter,
    Acts::CylinderVolumeBounds::Face barrelMaterialFace) {
  const auto assemblyElement = builder.findDetElementByName(assembly);
  if (!assemblyElement.has_value()) {
    throw std::runtime_error(
        std::format("Could not find assembly '{}'", assembly));
  }

  builder.barrelEndcap()
      .setAssembly(*assemblyElement)
      .setSensorAxes("XYZ", "XZY")
      .setLayerFilter(layerFilter)
      .onLayer(makeLayerCustomizer(builder, std::move(det), layerFilter,
                                   barrelMaterialFace))
      .onContainer([&builder](const dd4hep::DetElement& elem,
                              Acts::detail::ContainerNodePtr node)
                       -> Acts::detail::BlueprintNodePtr {
        node->setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
        node->setResizeStrategies(Acts::VolumeResizeStrategy::Gap,
                                  Acts::VolumeResizeStrategy::Gap);

        // This callback fires for every container node the barrelEndcap()
        // builder creates: each endcap sub-container, the barrel
        // sub-container, and the combined top-level Z-stack. Only the
        // *barrel* sub-container (e.g. "PixelBarrel") carries
        // negative/positive boundary material per the ODD convention; the
        // combined top-level node must be left alone, since its own
        // negative/positive discs get fused with its radial neighbor's cap
        // in the outer AxisR stack, and Acts refuses to fuse two portals
        // that both carry material.
        if (!std::string{builder.backend().nameOf(elem)}.ends_with("Barrel")) {
          return node;
        }

        // Container-level material: negative/positive disc faces only. The
        // outer face is deliberately left alone here -- it collides with the
        // outermost layer's own material designation on the fused portal
        // where this container meets its radial neighbor (see the Kategorie
        // 2 investigation).
        using enum Acts::CylinderVolumeBounds::Face;
        auto mat = std::make_shared<Acts::MaterialDesignatorBlueprintNode>(
            node->name() + "_boundary_mat");
        mat->configureFace(NegativeDisc,
                           Acts::AxisSpec::DeferredEquidistant(
                               kMatRBins, Acts::AxisDirection::AxisR),
                           Acts::AxisSpec::DeferredEquidistant(
                               kMatPhiBins, Acts::AxisDirection::AxisPhi));
        mat->configureFace(PositiveDisc,
                           Acts::AxisSpec::DeferredEquidistant(
                               kMatRBins, Acts::AxisDirection::AxisR),
                           Acts::AxisSpec::DeferredEquidistant(
                               kMatPhiBins, Acts::AxisDirection::AxisPhi));
        mat->addChild(std::move(node));
        return mat;
      })
      .addTo(outer);
}

void addDirectLayerGroupedSubsystem(
    const BlueprintBuilder& builder, Acts::ContainerBlueprintNode& outer,
    std::string assembly, std::string det, const std::regex& layerFilter,
    Acts::CylinderVolumeBounds::Face barrelMaterialFace) {
  const auto assemblyElement = builder.findDetElementByName(assembly);
  if (!assemblyElement.has_value()) {
    throw std::runtime_error(
        std::format("Could not find assembly '{}'", assembly));
  }

  auto barrels = builder.findBarrelElements(*assemblyElement);
  auto endcaps = builder.findEndcapElements(*assemblyElement);

  const std::string assemblyName{builder.backend().nameOf(*assemblyElement)};
  auto containerNode = std::make_shared<Acts::CylinderContainerBlueprintNode>(
      assemblyName, Acts::AxisDirection::AxisZ);

  auto layerCustomizer = makeLayerCustomizer(builder, std::move(det),
                                             layerFilter, barrelMaterialFace);

  auto sensorToLayerKey = [&](const dd4hep::DetElement& elem) {
    auto current = elem;
    const auto world = builder.backend().world();
    while (!(current == world)) {
      std::cmatch match;
      if (const std::string name{builder.backend().nameOf(current)};
          std::regex_search(name.c_str(), match, layerFilter) &&
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
    barrelNode->setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
    barrelNode->setResizeStrategies(Acts::VolumeResizeStrategy::Gap,
                                    Acts::VolumeResizeStrategy::Gap);
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
    endcapNode->setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
    endcapNode->setResizeStrategies(Acts::VolumeResizeStrategy::Gap,
                                    Acts::VolumeResizeStrategy::Gap);
    containerNode->addChild(std::move(endcapNode));
  }

  outer.addChild(std::move(containerNode));
}

}  // namespace

std::unique_ptr<Acts::TrackingGeometry> buildOpenDataDetectorBarrelEndcap(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
    const Acts::Logger& logger) {
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

  outer.addMaterial(
      "Beampipe_mat", [&](Acts::MaterialDesignatorBlueprintNode& mat) {
        using enum Acts::CylinderVolumeBounds::Face;
        mat.configureFace(
            OuterCylinder,
            Acts::AxisSpec::DeferredEquidistant(kMatPhiBins, AxisRPhi),
            Acts::AxisSpec::DeferredEquidistant(kMatZBins, AxisZ));
        mat.addChild(builder.backend().makeBeampipe());
      });

  // Per-subsystem barrel material face: matches the ODDs own convention
  // (Pixel layers carry material on their outer face, ShortStrips/LongStrips
  // on their inner face) so that no two radially-stacked layers independently
  // claim the same fused portal. Hardcoded here rather than read from the
  // DD4hep XML, since that annotation is ODD-specific and not guaranteed to
  // be available for other detector geometries.
  using enum Acts::CylinderVolumeBounds::Face;
  addBarrelEndcapSubsystem(builder, outer, "Pixels", "pix",
                           ActsPlugins::DD4hep::detail::kPixelLayerFilter,
                           OuterCylinder);
  addBarrelEndcapSubsystem(builder, outer, "ShortStrips", "ss",
                           ActsPlugins::DD4hep::detail::kShortStripLayerFilter,
                           InnerCylinder);
  addBarrelEndcapSubsystem(builder, outer, "LongStrips", "ls",
                           ActsPlugins::DD4hep::detail::kLongStripLayerFilter,
                           InnerCylinder);

  return root.construct(BlueprintOptions{}, gctx, logger);
}

std::unique_ptr<Acts::TrackingGeometry> buildOpenDataDetectorDirectLayer(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
    const Acts::Logger& logger) {
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

  outer.addMaterial(
      "Beampipe_mat", [&](Acts::MaterialDesignatorBlueprintNode& mat) {
        using enum Acts::CylinderVolumeBounds::Face;
        mat.configureFace(
            OuterCylinder,
            Acts::AxisSpec::DeferredEquidistant(kMatPhiBins, AxisRPhi),
            Acts::AxisSpec::DeferredEquidistant(kMatZBins, AxisZ));
        mat.addChild(builder.backend().makeBeampipe());
      });

  using enum Acts::CylinderVolumeBounds::Face;
  addDirectLayerSubsystem(builder, outer, "Pixels", "pix",
                          ActsPlugins::DD4hep::detail::kPixelLayerFilter,
                          OuterCylinder);
  addDirectLayerSubsystem(builder, outer, "ShortStrips", "ss",
                          ActsPlugins::DD4hep::detail::kShortStripLayerFilter,
                          InnerCylinder);
  addDirectLayerSubsystem(builder, outer, "LongStrips", "ls",
                          ActsPlugins::DD4hep::detail::kLongStripLayerFilter,
                          InnerCylinder);

  return root.construct(BlueprintOptions{}, gctx, logger);
}

std::unique_ptr<Acts::TrackingGeometry> buildOpenDataDetectorDirectLayerGrouped(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
    const Acts::Logger& logger) {
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

  outer.addMaterial(
      "Beampipe_mat", [&](Acts::MaterialDesignatorBlueprintNode& mat) {
        using enum Acts::CylinderVolumeBounds::Face;
        mat.configureFace(
            OuterCylinder,
            Acts::AxisSpec::DeferredEquidistant(kMatPhiBins, AxisRPhi),
            Acts::AxisSpec::DeferredEquidistant(kMatZBins, AxisZ));
        mat.addChild(builder.backend().makeBeampipe());
      });

  using enum Acts::CylinderVolumeBounds::Face;
  addDirectLayerGroupedSubsystem(builder, outer, "Pixels", "pix",
                                 ActsPlugins::DD4hep::detail::kPixelLayerFilter,
                                 OuterCylinder);
  addDirectLayerGroupedSubsystem(
      builder, outer, "ShortStrips", "ss",
      ActsPlugins::DD4hep::detail::kShortStripLayerFilter, InnerCylinder);
  addDirectLayerGroupedSubsystem(
      builder, outer, "LongStrips", "ls",
      ActsPlugins::DD4hep::detail::kLongStripLayerFilter, InnerCylinder);

  return root.construct(BlueprintOptions{}, gctx, logger);
}

}  // namespace ActsPlugins::DD4hep
