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

// Beampipe material, identical across all three construction methods.
void addBeampipe(const BlueprintBuilder& builder,
                 Acts::ContainerBlueprintNode& outer) {
  using enum Acts::CylinderVolumeBounds::Face;
  outer.addMaterial(
      "Beampipe_mat", [&](Acts::MaterialDesignatorBlueprintNode& mat) {
        mat.configureFace(OuterCylinder,
                          Acts::AxisSpec::DeferredEquidistant(
                              kMatPhiBins, Acts::AxisDirection::AxisRPhi),
                          Acts::AxisSpec::DeferredEquidistant(
                              kMatZBins, Acts::AxisDirection::AxisZ));
        mat.addChild(builder.backend().makeBeampipe());
      });
}

// A single standalone passive tube-shaped element (e.g. PST, Solenoid):
// looked up by name, materialized on the given faces, and added as a child
// of `outer`. A no-op if `elementName` does not exist in this geometry.
// Identical across all three construction methods.
void addPassiveCylinder(
    const BlueprintBuilder& builder, Acts::ContainerBlueprintNode& outer,
    const std::string& elementName,
    std::initializer_list<Acts::CylinderVolumeBounds::Face> faces) {
  const auto element = builder.findDetElementByName(elementName);
  if (!element.has_value()) {
    return;
  }
  outer.addMaterial(
      elementName + "_mat", [&](Acts::MaterialDesignatorBlueprintNode& mat) {
        for (const auto face : faces) {
          mat.configureFace(face,
                            Acts::AxisSpec::DeferredEquidistant(
                                kMatPhiBins, Acts::AxisDirection::AxisRPhi),
                            Acts::AxisSpec::DeferredEquidistant(
                                kMatZBins, Acts::AxisDirection::AxisZ));
        }
        mat.addChild(builder.backend().makePassiveCylinder(*element));
      });
}

// Pixel endplates: passive carbon-fiber discs beyond the outermost endcap
// layer on each side (z=-1950mm / +1950mm). In the ODD XML they share the
// same layer_pattern as the sensitive endcap layers
// (`PixelEndcapN\d|PixelEndplate`), but carry no sensitive modules of their
// own, so the generic sensor-layer discovery that populates `endcapNode`
// never finds them. Added here as an extra static child of whichever endcap
// element actually has one as a direct DD4hep child -- only
// PixelEndcapN/PixelEndcapP do, so this is a no-op everywhere else. Shared
// across all three construction methods, each of which builds its own
// per-endcap node differently but always has both the DD4hep endcap element
// and its resulting node available at the same point.
void addPixelEndplateIfPresent(const BlueprintBuilder& builder,
                               const dd4hep::DetElement& endcapElement,
                               Acts::BlueprintNode& endcapNode) {
  using enum Acts::CylinderVolumeBounds::Face;
  for (const auto& child : builder.backend().children(endcapElement)) {
    if (builder.backend().nameOf(child) != "PixelEndplate") {
      continue;
    }
    auto endplateMat = std::make_shared<Acts::MaterialDesignatorBlueprintNode>(
        endcapNode.name() + "_endplate_mat");
    endplateMat->configureFace(NegativeDisc,
                               Acts::AxisSpec::DeferredEquidistant(
                                   kMatRBins, Acts::AxisDirection::AxisR),
                               Acts::AxisSpec::DeferredEquidistant(
                                   kMatPhiBins, Acts::AxisDirection::AxisPhi));
    endplateMat->configureFace(PositiveDisc,
                               Acts::AxisSpec::DeferredEquidistant(
                                   kMatRBins, Acts::AxisDirection::AxisR),
                               Acts::AxisSpec::DeferredEquidistant(
                                   kMatPhiBins, Acts::AxisDirection::AxisPhi));
    endplateMat->addChild(builder.backend().makePassiveDisc(
        child, endcapNode.name() + "_PixelEndplate"));
    endcapNode.addChild(std::move(endplateMat));
  }
}

// Barrel container boundary_material: negative/positive disc faces on the
// "*Barrel" container as a whole, on top of (and independent from) whatever
// material its individual layers already carry. Per the ODD convention this
// only ever applies to barrel containers, never to endcaps -- their own
// Z-facing boundary is deliberately left unmaterialized at the container
// level, since it would otherwise collide with this same designation once
// barrel and endcaps are Z-stacked together and their touching faces get
// fused into one portal. Shared across all three construction methods.
Acts::detail::BlueprintNodePtr addBarrelBoundaryMaterial(
    Acts::detail::ContainerNodePtr node) {
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
}

void addDirectLayerSubsystem(
    const BlueprintBuilder& builder, Acts::ContainerBlueprintNode& outer,
    std::string assembly, std::string det, const std::regex& layerFilter,
    Acts::CylinderVolumeBounds::Face barrelMaterialFace,
    bool hasOuterBoundaryMaterial) {
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

  auto addLayerChildren = [&](const auto& elements, auto makeNode,
                              bool isBarrel) {
    for (const auto& element : elements) {
      auto node = makeNode(element);
      node->setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
      node->setResizeStrategies(Acts::VolumeResizeStrategy::Gap,
                                Acts::VolumeResizeStrategy::Gap);
      addPixelEndplateIfPresent(builder, element, *node);
      if (isBarrel) {
        containerNode->addChild(addBarrelBoundaryMaterial(std::move(node)));
      } else {
        containerNode->addChild(std::move(node));
      }
    }
  };

  addLayerChildren(
      barrels,
      [&](const auto& barrel) {
        return builder.layers()
            .barrel()
            .setSensorAxes("XYZ")
            .setLayerFilter(layerFilter)
            .setContainer(barrel)
            .onLayer(layerCustomizer)
            .build();
      },
      /*isBarrel=*/true);

  addLayerChildren(
      endcaps,
      [&](const auto& endcap) {
        return builder.layers()
            .endcap()
            .setSensorAxes("XZY")
            .setLayerFilter(layerFilter)
            .setContainer(endcap)
            .onLayer(layerCustomizer)
            .build();
      },
      /*isBarrel=*/false);

  if (!hasOuterBoundaryMaterial) {
    outer.addChild(std::move(containerNode));
    return;
  }

  using enum Acts::CylinderVolumeBounds::Face;
  auto outerMat = std::make_shared<Acts::MaterialDesignatorBlueprintNode>(
      assembly + "_outer_mat");
  outerMat->configureFace(OuterCylinder,
                          Acts::AxisSpec::DeferredEquidistant(
                              kMatPhiBins, Acts::AxisDirection::AxisRPhi),
                          Acts::AxisSpec::DeferredEquidistant(
                              kMatZBins, Acts::AxisDirection::AxisZ));
  outerMat->addChild(std::move(containerNode));
  outer.addChild(std::move(outerMat));
}

void addBarrelEndcapSubsystem(
    const BlueprintBuilder& builder, Acts::ContainerBlueprintNode& outer,
    std::string assembly, std::string det, const std::regex& layerFilter,
    Acts::CylinderVolumeBounds::Face barrelMaterialFace,
    bool hasOuterBoundaryMaterial) {
  const auto assemblyElement = builder.findDetElementByName(assembly);
  if (!assemblyElement.has_value()) {
    throw std::runtime_error(
        std::format("Could not find assembly '{}'", assembly));
  }

  auto topNode =
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

            addPixelEndplateIfPresent(builder, elem, *node);

            // This callback fires for every container node the barrelEndcap()
            // builder creates: each endcap sub-container, the barrel
            // sub-container, and the combined top-level Z-stack. Only the
            // *barrel* sub-container (e.g. "PixelBarrel") carries
            // negative/positive boundary material per the ODD convention; the
            // combined top-level node must be left alone, since its own
            // negative/positive discs get fused with its radial neighbor's cap
            // in the outer AxisR stack, and Acts refuses to fuse two portals
            // that both carry material. The outer cylinder face is handled
            // separately below, on the fully Z-stacked top-level node returned
            // by build() -- not here, since materializing it on this "Barrel"
            // sub-container (before it gets Z-stacked together with the
            // endcaps) collides with the subsequent radius-unification across
            // the whole Z-stack (see the Kategorie 2 investigation).
            if (!std::string{builder.backend().nameOf(elem)}.ends_with(
                    "Barrel")) {
              return node;
            }
            return addBarrelBoundaryMaterial(std::move(node));
          })
          .build();

  if (!hasOuterBoundaryMaterial) {
    outer.addChild(std::move(topNode));
    return;
  }

  // Outer boundary_material: wrap the fully Z-stacked, radius-unified
  // top-level node returned by build() -- not accessible via onContainer,
  // since BarrelEndcapAssembler constructs this node directly and adds it as
  // a child without ever routing it back through that callback. Doing it
  // here, after the Z-stack's radius has already been unified across
  // barrel+endcaps, avoids the "merged when stacking" conflict that occurs
  // when materializing the outer face on a raw pre-stack constituent (see
  // the Kategorie 2 investigation).
  using enum Acts::CylinderVolumeBounds::Face;
  auto outerMat = std::make_shared<Acts::MaterialDesignatorBlueprintNode>(
      assembly + "_outer_mat");
  outerMat->configureFace(OuterCylinder,
                          Acts::AxisSpec::DeferredEquidistant(
                              kMatPhiBins, Acts::AxisDirection::AxisRPhi),
                          Acts::AxisSpec::DeferredEquidistant(
                              kMatZBins, Acts::AxisDirection::AxisZ));
  outerMat->addChild(std::move(topNode));
  outer.addChild(std::move(outerMat));
}

void addDirectLayerGroupedSubsystem(
    const BlueprintBuilder& builder, Acts::ContainerBlueprintNode& outer,
    std::string assembly, std::string det, const std::regex& layerFilter,
    Acts::CylinderVolumeBounds::Face barrelMaterialFace,
    bool hasOuterBoundaryMaterial) {
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
    containerNode->addChild(addBarrelBoundaryMaterial(std::move(barrelNode)));
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
    addPixelEndplateIfPresent(builder, endcap, *endcapNode);
    containerNode->addChild(std::move(endcapNode));
  }

  if (!hasOuterBoundaryMaterial) {
    outer.addChild(std::move(containerNode));
    return;
  }

  using enum Acts::CylinderVolumeBounds::Face;
  auto outerMat = std::make_shared<Acts::MaterialDesignatorBlueprintNode>(
      assembly + "_outer_mat");
  outerMat->configureFace(OuterCylinder,
                          Acts::AxisSpec::DeferredEquidistant(
                              kMatPhiBins, Acts::AxisDirection::AxisRPhi),
                          Acts::AxisSpec::DeferredEquidistant(
                              kMatZBins, Acts::AxisDirection::AxisZ));
  outerMat->addChild(std::move(containerNode));
  outer.addChild(std::move(outerMat));
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

  addBeampipe(builder, outer);

  // Per-subsystem barrel material face: matches the ODDs own convention
  // (Pixel layers carry material on their outer face, ShortStrips/LongStrips
  // on their inner face) so that no two radially-stacked layers independently
  // claim the same fused portal. Hardcoded here rather than read from the
  // DD4hep XML, since that annotation is ODD-specific and not guaranteed to
  // be available for other detector geometries.
  using enum Acts::CylinderVolumeBounds::Face;
  // PixelBarrel's own "outer" boundary_material flag shares the exact same
  // binning constants (mat_pix_barrel_bPhi/bZ) as every individual Pixel
  // layer's own "outer" layer_material flag -- unlike the LongStrips/Solenoid
  // case, this isn't two distinct designations that happen to coincide, it's
  // the same one, re-stated at the container level. Confirmed in Gen1: no
  // separate container-boundary surface exists anywhere near the Pixel
  // barrel's outer radius, only the outermost layer's own material. Enabling
  // it here would double-count that layer's material at a second, nearby but
  // distinct radius instead of reproducing Gen1's single merged surface.
  addBarrelEndcapSubsystem(builder, outer, "Pixels", "pix",
                           ActsPlugins::DD4hep::detail::kPixelLayerFilter,
                           OuterCylinder, /*hasOuterBoundaryMaterial=*/false);

  // Passive Support Tube (PST): a thin carbon-fiber support cylinder between
  // the Pixel and ShortStrips subsystems. It is not a tracker sub-detector,
  // just a single passive tube-shaped element with its own material.
  addPassiveCylinder(builder, outer, "PST", {OuterCylinder});

  addBarrelEndcapSubsystem(builder, outer, "ShortStrips", "ss",
                           ActsPlugins::DD4hep::detail::kShortStripLayerFilter,
                           InnerCylinder, /*hasOuterBoundaryMaterial=*/true);
  // LongStripBarrel has no "outer" boundary_material flag of its own in the
  // XML, but that is only because the adjacent Solenoid's own "inner" flag
  // is meant to cover this shared boundary instead -- an assumption that
  // holds in Gen1, where the two volumes are snapped together with no gap,
  // but not in Gen3, where VolumeAttachmentStrategy::Gap keeps them apart at
  // their true, distinct radii (see the LongStrips/Solenoid investigation).
  // Materializing it here restores full coverage of that boundary.
  addBarrelEndcapSubsystem(builder, outer, "LongStrips", "ls",
                           ActsPlugins::DD4hep::detail::kLongStripLayerFilter,
                           InnerCylinder, /*hasOuterBoundaryMaterial=*/true);

  // Solenoid: a passive aluminum tube-shaped element outside LongStrips,
  // structurally identical to the beampipe/PST case. Unlike PST, its ODD XML
  // carries both a `layer_material surface="representing"` (-> outer face,
  // same convention as PST/beampipe) and a `boundary_material surface="inner"`
  // (-> inner face) on the same element.
  addPassiveCylinder(builder, outer, "Solenoid",
                     {InnerCylinder, OuterCylinder});

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

  addBeampipe(builder, outer);

  using enum Acts::CylinderVolumeBounds::Face;
  addDirectLayerSubsystem(builder, outer, "Pixels", "pix",
                          ActsPlugins::DD4hep::detail::kPixelLayerFilter,
                          OuterCylinder, /*hasOuterBoundaryMaterial=*/false);
  addPassiveCylinder(builder, outer, "PST", {OuterCylinder});
  addDirectLayerSubsystem(builder, outer, "ShortStrips", "ss",
                          ActsPlugins::DD4hep::detail::kShortStripLayerFilter,
                          InnerCylinder, /*hasOuterBoundaryMaterial=*/true);
  addDirectLayerSubsystem(builder, outer, "LongStrips", "ls",
                          ActsPlugins::DD4hep::detail::kLongStripLayerFilter,
                          InnerCylinder, /*hasOuterBoundaryMaterial=*/true);
  addPassiveCylinder(builder, outer, "Solenoid",
                     {InnerCylinder, OuterCylinder});

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

  addBeampipe(builder, outer);

  using enum Acts::CylinderVolumeBounds::Face;
  addDirectLayerGroupedSubsystem(builder, outer, "Pixels", "pix",
                                 ActsPlugins::DD4hep::detail::kPixelLayerFilter,
                                 OuterCylinder,
                                 /*hasOuterBoundaryMaterial=*/false);
  addPassiveCylinder(builder, outer, "PST", {OuterCylinder});
  addDirectLayerGroupedSubsystem(
      builder, outer, "ShortStrips", "ss",
      ActsPlugins::DD4hep::detail::kShortStripLayerFilter, InnerCylinder,
      /*hasOuterBoundaryMaterial=*/true);
  addDirectLayerGroupedSubsystem(
      builder, outer, "LongStrips", "ls",
      ActsPlugins::DD4hep::detail::kLongStripLayerFilter, InnerCylinder,
      /*hasOuterBoundaryMaterial=*/true);
  addPassiveCylinder(builder, outer, "Solenoid",
                     {InnerCylinder, OuterCylinder});

  return root.construct(BlueprintOptions{}, gctx, logger);
}

}  // namespace ActsPlugins::DD4hep
