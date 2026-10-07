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

using Face = Acts::CylinderVolumeBounds::Face;

// Material bin counts, chosen analogous to Gen1. Phi is 144 everywhere in Gen1,
// so one shared constant covers it; Z/R differ per subsystem.
constexpr std::size_t kMatPhiBins = 144;

// beam pipe
constexpr std::size_t kBeampipeZBins = 10;

// pixel
constexpr std::size_t kPixelBarrelZBins = 200;
constexpr std::size_t kPixelEndcapRBins = 150;

// PST
constexpr std::size_t kPstZBins = 100;

// short strip
constexpr std::size_t kShortStripBarrelZBins = 150;
constexpr std::size_t kShortStripEndcapRBins = 150;

// long strip
constexpr std::size_t kLongStripBarrelZBins = 250;
constexpr std::size_t kLongStripEndcapRBins = 100;

// Solenoid
constexpr std::size_t kSolenoidZBins = 100;

// Shared by the Pixel/ShortStrip/LongStrip "*Barrel" containers' own
// negative/positive boundary_material in Gen1.
constexpr std::size_t kBarrelContainerBoundaryRBins = 150;

// LongStrips has no "outer" boundary_material of its own in Gen1 -- the
// Solenoid's "inner" flag covers it instead, so Gen3's equivalent addition
// reuses the Solenoid's bin count.
constexpr std::size_t kLongStripOuterBoundaryZBins = kSolenoidZBins;

// Calorimeter collector (no Gen1 equivalent). It carries nearly all of the
// calorimeter material, which changes steeply across the barrel-endcap
// transition (1 < |eta| < 2), so it needs fine bins: 100 z bins over the
// resized collector length give ~80mm bins, i.e. < 0.03 in eta at |eta| ~ 1.5.
constexpr std::size_t kCaloCollectorZBins = 100;
// The disc collectors carry the endcap calorimeter material. Its inner edge
// (0 -> ~105 X0 between r ~205mm and ~320mm, i.e. 3.0 < |eta| < 3.45) needs
// fine radial bins: 400 bins over the disc give ~3mm. The calorimeter is
// nearly phi-symmetric, so 36 phi bins suffice and keep the bin count (and
// thus the statistics per bin) at that of 144 x 100.
constexpr std::size_t kCaloCollectorDiscRBins = 400;
constexpr std::size_t kCaloCollectorDiscPhiBins = 36;

// Short, stable label for a face, used to build a per-face material key (a
// designator node can configure more than one face, e.g. the Solenoid's
// Inner+Outer cylinder, or an endcap layer's Negative+PositiveDisc -- the
// node's own name alone would not be unique in that case).
std::string faceLabel(Face face) {
  switch (face) {
    case Face::NegativeDisc:
      return "NegativeDisc";
    case Face::PositiveDisc:
      return "PositiveDisc";
    case Face::OuterCylinder:
      return "OuterCylinder";
    case Face::InnerCylinder:
      return "InnerCylinder";
    default:
      throw std::invalid_argument("faceLabel: unsupported face");
  }
}

// Configures `face` as a cylinder mantle: bins in (RPhi, Z). Used wherever a
// thin cylindrical shell carries material on one of its two mantle faces --
// the beampipe/PST/Solenoid, barrel layers, and outer subsystem boundaries.
// Stable key: derived from the designator node's own name plus the face, so
// a map stays valid across GeometryIdentifier churn from geometry changes.
void configureCylinderFace(Acts::MaterialDesignatorBlueprintNode& mat,
                           Face face, std::size_t phiBins, std::size_t zBins) {
  mat.configureFace(
      face,
      Acts::AxisSpec::DeferredEquidistant(phiBins,
                                          Acts::AxisDirection::AxisRPhi),
      Acts::AxisSpec::DeferredEquidistant(zBins, Acts::AxisDirection::AxisZ),
      mat.name() + "_" + faceLabel(face));
}

// Configures `face` as a flat disc: bins in (R, Phi). Used for endcap layer
// material, pixel endplates, and container-level Negative/PositiveDisc
// boundary material. Stable key: see @ref configureCylinderFace.
void configureDiscFace(Acts::MaterialDesignatorBlueprintNode& mat, Face face,
                       std::size_t phiBins, std::size_t rBins) {
  mat.configureFace(
      face,
      Acts::AxisSpec::DeferredEquidistant(rBins, Acts::AxisDirection::AxisR),
      Acts::AxisSpec::DeferredEquidistant(phiBins,
                                          Acts::AxisDirection::AxisPhi),
      mat.name() + "_" + faceLabel(face));
}

// Every subsystem container in this file uses the same Gap attachment/resize
// strategy, so that radially or longitudinally adjacent volumes are kept
// apart at their true nominal extents (via a genuine gap-filler volume)
// instead of being snapped together with zero clearance.
void setGapAttachment(const Acts::detail::ContainerNodePtr& node) {
  node->setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
  node->setResizeStrategies(Acts::VolumeResizeStrategy::Gap,
                            Acts::VolumeResizeStrategy::Gap);
}

// Looks up `assembly` by name, throwing if this geometry has none.
dd4hep::DetElement findAssemblyOrThrow(const BlueprintBuilder& builder,
                                       const std::string& assembly) {
  const auto assemblyElement = builder.findDetElementByName(assembly);
  if (!assemblyElement.has_value()) {
    throw std::runtime_error(
        std::format("Could not find assembly '{}'", assembly));
  }
  return *assemblyElement;
}

auto makeLayerCustomizer(const BlueprintBuilder& builder, std::string det,
                         std::regex layerFilter, Face barrelMaterialFace,
                         std::size_t barrelZBins, std::size_t endcapRBins) {
  return [&builder, det = std::move(det), layerFilter = std::move(layerFilter),
          barrelMaterialFace, barrelZBins,
          endcapRBins](const std::optional<dd4hep::DetElement>& elem,
                       Acts::detail::LayerNodePtr layer)
             -> Acts::detail::BlueprintNodePtr {
    layer->setEnvelope(detail::kLayerEnvelope);

    const std::string elemName =
        elem.has_value() ? builder.backend().nameOf(*elem) : layer->name();
    const int layerIdx = detail::layerIndexFromName(elemName, layerFilter);

    using SrfArrayNavPol = Acts::SurfaceArrayNavigationPolicy;
    using enum SrfArrayNavPol::LayerType;

    SrfArrayNavPol::Config navCfg;
    navCfg.envelope = detail::kLayerEnvelope;

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
      configureCylinderFace(*matNode, barrelMaterialFace, kMatPhiBins,
                            barrelZBins);
    } else {
      // Endcap layer
      navCfg.layerType = Disc;
      navCfg.bins = {builder.backend().constant("{}_e_sf_b_r", det),
                     builder.backend().constant("{}_e_sf_b_phi", det)};

      // Endcap: thin disc "pancake" -> both flat faces carry material
      configureDiscFace(*matNode, Face::NegativeDisc, kMatPhiBins, endcapRBins);
      configureDiscFace(*matNode, Face::PositiveDisc, kMatPhiBins, endcapRBins);
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
  outer.addMaterial("Beampipe_mat",
                    [&](Acts::MaterialDesignatorBlueprintNode& mat) {
                      configureCylinderFace(mat, Face::OuterCylinder,
                                            kMatPhiBins, kBeampipeZBins);
                      mat.addChild(builder.backend().makeBeampipe());
                    });
}

// A single standalone passive tube-shaped element (e.g. PST, Solenoid):
// looked up by name, materialized on the given faces, and added as a child
// of `outer`. A no-op if `elementName` does not exist in this geometry.
// Identical across all three construction methods.
void addPassiveCylinder(const BlueprintBuilder& builder,
                        Acts::ContainerBlueprintNode& outer,
                        const std::string& elementName,
                        std::initializer_list<Face> faces, std::size_t zBins) {
  const auto element = builder.findDetElementByName(elementName);
  if (!element.has_value()) {
    return;
  }
  outer.addMaterial(
      elementName + "_mat", [&](Acts::MaterialDesignatorBlueprintNode& mat) {
        for (const auto face : faces) {
          configureCylinderFace(mat, face, kMatPhiBins, zBins);
        }
        mat.addChild(builder.backend().makePassiveElement(*element));
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
  for (const auto& child : builder.backend().children(endcapElement)) {
    if (builder.backend().nameOf(child) != "PixelEndplate") {
      continue;
    }
    auto endplateMat = std::make_shared<Acts::MaterialDesignatorBlueprintNode>(
        endcapNode.name() + "_endplate_mat");
    configureDiscFace(*endplateMat, Face::NegativeDisc, kMatPhiBins,
                      kPixelEndcapRBins);
    configureDiscFace(*endplateMat, Face::PositiveDisc, kMatPhiBins,
                      kPixelEndcapRBins);
    endplateMat->addChild(builder.backend().makePassiveElement(
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
  auto mat = std::make_shared<Acts::MaterialDesignatorBlueprintNode>(
      node->name() + "_boundary_mat");
  configureDiscFace(*mat, Face::NegativeDisc, kMatPhiBins,
                    kBarrelContainerBoundaryRBins);
  configureDiscFace(*mat, Face::PositiveDisc, kMatPhiBins,
                    kBarrelContainerBoundaryRBins);
  mat->addChild(std::move(node));
  return mat;
}

// Outer boundary_material: the OuterCylinder face of a fully assembled
// subsystem, added to `outer`. Must run on the fully Z-stacked, radius-
// unified top-level node -- never on one of its still-unstacked constituents.
// A constituent's outer face is later merged with its Z-neighbors' into the
// larger unified surface, and Blueprint construction aborts when a portal
// face carrying designated material has to be merged.
// For BarrelEndcap this cannot happen in the onContainer callback either:
// BarrelEndcapAssembler creates the top-level node itself and never routes it
// through that callback, so the material is attached to the node returned by
// build(). Shared across all three construction methods.
void addOuterBoundaryMaterial(Acts::ContainerBlueprintNode& outer,
                              Acts::detail::BlueprintNodePtr node,
                              const std::string& assembly, std::size_t zBins) {
  outer.addMaterial(
      assembly + "_outer_mat", [&](Acts::MaterialDesignatorBlueprintNode& mat) {
        configureCylinderFace(mat, Face::OuterCylinder, kMatPhiBins, zBins);
        mat.addChild(std::move(node));
      });
}

// Synthetic "calorimeter material collector" barrel: a thin annular cylinder
// with no corresponding DD4hep element, sitting in the radial gap between
// the Solenoid (rmax=1200mm) and the ECal barrel (rmin=1250mm) -- only 50mm
// wide in total, so the collector is centered in it with as much clearance
// as that gap allows on either side (~15mm to the Solenoid, ~25mm to the
// ECal). Added as one more radial shell of the tracker stack, exactly like
// @ref addPassiveCylinder adds the Solenoid. Its z half-length (3154mm) is
// deliberately the longest of all tracker shells, so that it sets the
// tracker's z-extent: it lies beyond every tracker envelope (|z| <= 3150mm),
// so every subsystem is lengthened with gap volumes at its ends (Gap resize
// strategy). Otherwise the longest subsystem's end faces would be the
// material faces of its outermost endcap layer, which can't be merged with
// the other shells' end faces during the radial stacking. The matching two
// disc collectors sit just beyond, see @ref addCaloMaterialCollectorDisc.
// The exact value 1215mm * 3180mm / 1225mm = 3154mm makes the line from the
// origin through the end of the material face (r=1215mm) hit the disc
// collectors' material faces (|z|=3180mm) exactly at their outer edge
// (r=1225mm): tracks from the origin then cross exactly one of the two
// collectors, with neither overlap (double counting) nor gap.
// Material is designated on the inner face (r=1215mm): the outer face is
// the outer face of the whole tracker, which is merged with the disc
// collectors' outer faces when they are stacked along z.
//
// Acts::MaterialInteractionAssignment::assign() attributes every recorded
// material interaction to its single nearest materialized surface. Without
// this collector, interactions in the tracker-calorimeter gap (cables,
// support structure, or calorimeter material reaching slightly inward) would
// be misattributed to the Solenoid's own outer face or another nearby
// tracker surface instead. Shared across all three construction methods.
void addCaloMaterialCollector(const BlueprintBuilder& builder,
                              Acts::ContainerBlueprintNode& outer) {
  constexpr double kRMin = 1215 * Acts::UnitConstants::mm;
  constexpr double kRMax = 1225 * Acts::UnitConstants::mm;
  constexpr double kHalfZ = 3154 * Acts::UnitConstants::mm;

  outer.addMaterial(
      "CaloMaterialCollectorBarrel_mat",
      [&](Acts::MaterialDesignatorBlueprintNode& mat) {
        configureCylinderFace(mat, Face::InnerCylinder, kMatPhiBins,
                              kCaloCollectorZBins);
        mat.addChild(builder.backend().makeMaterialCollector(
            kRMin, kRMax, kHalfZ, 0.0, "CaloMaterialCollectorBarrel"));
      });
}

// AxisR stack holding every tracker shell (Pixels ... Solenoid and the
// calorimeter collector), but not the beampipe. Kept separate from the
// beampipe so that its z-extent is set by the tracker itself, not by the
// much longer beampipe (see @ref finalizeOuter).
std::shared_ptr<Acts::CylinderContainerBlueprintNode> makeTrackerNode() {
  auto tracker = std::make_shared<Acts::CylinderContainerBlueprintNode>(
      "Tracker", Acts::AxisDirection::AxisR);
  tracker->setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
  return tracker;
}

// Endcap calorimeter material collector: a thin disc volume on each side of
// the tracker, in the gap between the end of the tracker envelopes
// (|z| = 3150mm) and the front of the ECal endcap (|z| = 3200mm). Its radial
// extent is nominal, the AxisZ stack in @ref finalizeOuter resizes it to the
// tracker's radial extent. Material is designated on the face pointing away
// from the tracker (towards the calorimeter), which is only fused with the
// adjacent gap volume and never merged.
void addCaloMaterialCollectorDisc(const BlueprintBuilder& builder,
                                  Acts::ContainerBlueprintNode& trackerZ,
                                  bool positive) {
  constexpr double kRMin = 30 * Acts::UnitConstants::mm;
  constexpr double kRMax = 1225 * Acts::UnitConstants::mm;
  constexpr double kHalfZ = 5 * Acts::UnitConstants::mm;
  constexpr double kAbsZ = 3175 * Acts::UnitConstants::mm;

  const std::string name =
      positive ? "CaloMaterialCollectorDiscP" : "CaloMaterialCollectorDiscN";
  trackerZ.addMaterial(
      name + "_mat", [&](Acts::MaterialDesignatorBlueprintNode& mat) {
        configureDiscFace(mat,
                          positive ? Face::PositiveDisc : Face::NegativeDisc,
                          kCaloCollectorDiscPhiBins, kCaloCollectorDiscRBins);
        mat.addChild(builder.backend().makeMaterialCollector(
            kRMin, kRMax, kHalfZ, positive ? kAbsZ : -kAbsZ, name));
      });
}

// Assembles the top level and returns the constructed geometry:
//
//   root
//   └─ outer (AxisR)
//      ├─ beampipe                         (already added to outer)
//      └─ TrackerZ (AxisZ, Gap attachment and resize)
//         ├─ collector disc, -z
//         ├─ tracker (AxisR)
//         └─ collector disc, +z
//
// In an AxisR stack every shell is resized to the common z-extent, which is
// set by the beampipe (+-4m). If the tracker shells were stacked with the
// beampipe directly, they would all be stretched to +-4m: the outer faces of
// ShortStrips/Solenoid would extend through the endcap calorimeter region and
// collector discs on the stack's end faces would sit behind the front of the
// ECal endcap (z=3.2m), so the endcap calorimeter material would be
// attributed to the nearer tracker endcap discs. Wrapping the tracker in an
// AxisZ container with Gap resize strategy fills the extra length with gap
// volumes instead, so the tracker keeps its own z-extent and the collector
// discs sit in front of the ECal endcap.
//
// The collector discs are separate volumes since material can't be designated
// on the end faces of the tracker stack itself: those faces are composed of
// its children's faces, which are merged during the radial stacking.
// Together with the collector cylinder the discs form a closed
// two-discs-plus-one-barrel cage around the tracker. Shared across all three
// construction methods.
std::unique_ptr<Acts::TrackingGeometry> finalizeOuter(
    const BlueprintBuilder& builder, Acts::Blueprint& root,
    std::shared_ptr<Acts::CylinderContainerBlueprintNode> outer,
    std::shared_ptr<Acts::CylinderContainerBlueprintNode> tracker,
    const Acts::GeometryContext& gctx, const Acts::Logger& logger) {
  auto trackerZ = std::make_shared<Acts::CylinderContainerBlueprintNode>(
      "TrackerZ", Acts::AxisDirection::AxisZ);
  trackerZ->setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
  trackerZ->setResizeStrategies(Acts::VolumeResizeStrategy::Gap,
                                Acts::VolumeResizeStrategy::Gap);
  addCaloMaterialCollectorDisc(builder, *trackerZ, /*positive=*/false);
  trackerZ->addChild(std::move(tracker));
  addCaloMaterialCollectorDisc(builder, *trackerZ, /*positive=*/true);

  outer->addChild(std::move(trackerZ));
  root.addChild(std::move(outer));
  return root.construct(Acts::BlueprintOptions{}, gctx, logger);
}

// Adds `node` (either a barrel or an endcap sub-container) as a Z-stacked
// child of `containerNode`, applying the uniform Gap attachment/resize
// strategy plus whichever additional per-role material a barrel or endcap
// container needs (see @ref addBarrelBoundaryMaterial and @ref
// addPixelEndplateIfPresent). Shared between the two flat-layer construction
// methods (DirectLayer / DirectLayerGrouped).
void addSubsystemChild(const BlueprintBuilder& builder,
                       const dd4hep::DetElement& element,
                       Acts::ContainerBlueprintNode& containerNode,
                       Acts::detail::ContainerNodePtr node, bool isBarrel) {
  setGapAttachment(node);
  if (isBarrel) {
    containerNode.addChild(addBarrelBoundaryMaterial(std::move(node)));
    return;
  }
  addPixelEndplateIfPresent(builder, element, *node);
  containerNode.addChild(std::move(node));
}

void addDirectLayerSubsystem(const BlueprintBuilder& builder,
                             Acts::ContainerBlueprintNode& outer,
                             std::string assembly, std::string det,
                             const std::regex& layerFilter,
                             Face barrelMaterialFace, std::size_t barrelZBins,
                             std::size_t endcapRBins,
                             std::optional<std::size_t> outerBoundaryZBins) {
  const auto assemblyElement = findAssemblyOrThrow(builder, assembly);
  auto barrels = builder.findBarrelElements(assemblyElement);
  auto endcaps = builder.findEndcapElements(assemblyElement);

  auto containerNode = std::make_shared<Acts::CylinderContainerBlueprintNode>(
      builder.backend().nameOf(assemblyElement), Acts::AxisDirection::AxisZ);

  auto layerCustomizer =
      makeLayerCustomizer(builder, std::move(det), layerFilter,
                          barrelMaterialFace, barrelZBins, endcapRBins);

  for (const auto& barrel : barrels) {
    auto node = builder.layers()
                    .barrel()
                    .setSensorAxes("XYZ")
                    .setLayerFilter(layerFilter)
                    .setContainer(barrel)
                    .onLayer(layerCustomizer)
                    .build();
    addSubsystemChild(builder, barrel, *containerNode, std::move(node),
                      /*isBarrel=*/true);
  }

  for (const auto& endcap : endcaps) {
    auto node = builder.layers()
                    .endcap()
                    .setSensorAxes("XZY")
                    .setLayerFilter(layerFilter)
                    .setContainer(endcap)
                    .onLayer(layerCustomizer)
                    .build();
    addSubsystemChild(builder, endcap, *containerNode, std::move(node),
                      /*isBarrel=*/false);
  }

  if (!outerBoundaryZBins.has_value()) {
    outer.addChild(std::move(containerNode));
    return;
  }
  addOuterBoundaryMaterial(outer, std::move(containerNode), assembly,
                           *outerBoundaryZBins);
}

void addBarrelEndcapSubsystem(const BlueprintBuilder& builder,
                              Acts::ContainerBlueprintNode& outer,
                              std::string assembly, std::string det,
                              const std::regex& layerFilter,
                              Face barrelMaterialFace, std::size_t barrelZBins,
                              std::size_t endcapRBins,
                              std::optional<std::size_t> outerBoundaryZBins) {
  const auto assemblyElement = findAssemblyOrThrow(builder, assembly);

  auto topNode =
      builder.barrelEndcap()
          .setAssembly(assemblyElement)
          .setSensorAxes("XYZ", "XZY")
          .setLayerFilter(layerFilter)
          .onLayer(makeLayerCustomizer(builder, std::move(det), layerFilter,
                                       barrelMaterialFace, barrelZBins,
                                       endcapRBins))
          .onContainer([&builder](const dd4hep::DetElement& elem,
                                  Acts::detail::ContainerNodePtr node)
                           -> Acts::detail::BlueprintNodePtr {
            setGapAttachment(node);
            addPixelEndplateIfPresent(builder, elem, *node);

            // This callback fires for every container node the barrelEndcap()
            // builder creates: each endcap sub-container, the barrel
            // sub-container, and the combined top-level Z-stack. Only the
            // *barrel* sub-container (e.g. "PixelBarrel") carries
            // negative/positive boundary material per the ODD convention; the
            // combined top-level node must be left alone, since its own
            // negative/positive discs get fused with its radial neighbor's
            // cap in the outer AxisR stack, and Acts refuses to fuse two
            // portals that both carry material. The outer cylinder face is
            // handled separately below, on the fully Z-stacked top-level node
            // returned by build() -- not here, since on this "Barrel"
            // sub-container the face would later be merged with the endcaps'
            // outer faces, which aborts construction for a face carrying
            // material (see @ref addOuterBoundaryMaterial).
            if (!builder.backend().nameOf(elem).ends_with("Barrel")) {
              return node;
            }
            return addBarrelBoundaryMaterial(std::move(node));
          })
          .build();

  if (!outerBoundaryZBins.has_value()) {
    outer.addChild(std::move(topNode));
    return;
  }
  addOuterBoundaryMaterial(outer, std::move(topNode), assembly,
                           *outerBoundaryZBins);
}

void addDirectLayerGroupedSubsystem(
    const BlueprintBuilder& builder, Acts::ContainerBlueprintNode& outer,
    std::string assembly, std::string det, const std::regex& layerFilter,
    Face barrelMaterialFace, std::size_t barrelZBins, std::size_t endcapRBins,
    std::optional<std::size_t> outerBoundaryZBins) {
  const auto assemblyElement = findAssemblyOrThrow(builder, assembly);
  auto barrels = builder.findBarrelElements(assemblyElement);
  auto endcaps = builder.findEndcapElements(assemblyElement);

  auto containerNode = std::make_shared<Acts::CylinderContainerBlueprintNode>(
      builder.backend().nameOf(assemblyElement), Acts::AxisDirection::AxisZ);

  auto layerCustomizer =
      makeLayerCustomizer(builder, std::move(det), layerFilter,
                          barrelMaterialFace, barrelZBins, endcapRBins);

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
    auto node = builder.layersFromSensors()
                    .barrel()
                    .setSensorAxes("XYZ")
                    .setSensors(std::move(sensors))
                    .setContainerName(builder.backend().nameOf(barrel))
                    .groupBy(sensorToLayerKey)
                    .onLayer(layerCustomizer)
                    .build();
    addSubsystemChild(builder, barrel, *containerNode, std::move(node),
                      /*isBarrel=*/true);
  }

  for (const auto& endcap : endcaps) {
    auto sensors = builder.resolveSensitives(endcap);
    auto node = builder.layersFromSensors()
                    .endcap()
                    .setSensorAxes("XZY")
                    .setSensors(std::move(sensors))
                    .setContainerName(builder.backend().nameOf(endcap))
                    .groupBy(sensorToLayerKey)
                    .onLayer(layerCustomizer)
                    .build();
    addSubsystemChild(builder, endcap, *containerNode, std::move(node),
                      /*isBarrel=*/false);
  }

  if (!outerBoundaryZBins.has_value()) {
    outer.addChild(std::move(containerNode));
    return;
  }
  addOuterBoundaryMaterial(outer, std::move(containerNode), assembly,
                           *outerBoundaryZBins);
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

  auto outerNode = std::make_shared<CylinderContainerBlueprintNode>(
      "OpenDataDetector", AxisR);
  outerNode->setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
  auto& outer = *outerNode;

  addBeampipe(builder, outer);

  auto trackerNode = makeTrackerNode();
  auto& tracker = *trackerNode;

  // Per-subsystem barrel material face: matches the ODDs own convention
  // (Pixel layers carry material on their outer face, ShortStrips/LongStrips
  // on their inner face) so that no two radially-stacked layers independently
  // claim the same fused portal. Hardcoded here rather than read from the
  // DD4hep XML, since that annotation is ODD-specific and not guaranteed to
  // be available for other detector geometries.
  using enum Face;
  // PixelBarrel's own "outer" boundary_material flag shares the exact same
  // binning constants (mat_pix_barrel_bPhi/bZ) as every individual Pixel
  // layer's own "outer" layer_material flag -- unlike the LongStrips/Solenoid
  // case, this isn't two distinct designations that happen to coincide, it's
  // the same one, re-stated at the container level. Confirmed in Gen1: no
  // separate container-boundary surface exists anywhere near the Pixel
  // barrel's outer radius, only the outermost layer's own material. Enabling
  // it here would double-count that layer's material at a second, nearby but
  // distinct radius instead of reproducing Gen1's single merged surface.
  addBarrelEndcapSubsystem(builder, tracker, "Pixels", "pix",
                           ActsPlugins::DD4hep::detail::kPixelLayerFilter,
                           OuterCylinder, kPixelBarrelZBins, kPixelEndcapRBins,
                           /*outerBoundaryZBins=*/std::nullopt);

  // Passive Support Tube (PST): a thin carbon-fiber support cylinder between
  // the Pixel and ShortStrips subsystems. It is not a tracker sub-detector,
  // just a single passive tube-shaped element with its own material.
  addPassiveCylinder(builder, tracker, "PST", {OuterCylinder}, kPstZBins);

  addBarrelEndcapSubsystem(builder, tracker, "ShortStrips", "ss",
                           ActsPlugins::DD4hep::detail::kShortStripLayerFilter,
                           InnerCylinder, kShortStripBarrelZBins,
                           kShortStripEndcapRBins,
                           /*outerBoundaryZBins=*/kShortStripBarrelZBins);
  // LongStripBarrel has no "outer" boundary_material flag of its own in the
  // XML, but that is only because the adjacent Solenoid's own "inner" flag
  // is meant to cover this shared boundary instead -- an assumption that
  // holds in Gen1, where the two volumes are snapped together with no gap,
  // but not in Gen3, where VolumeAttachmentStrategy::Gap keeps them apart at
  // their true, distinct radii (see the LongStrips/Solenoid investigation).
  // Materializing it here restores full coverage of that boundary, reusing
  // the Solenoid's own bin count since that's the designation being
  // reproduced (see kLongStripOuterBoundaryZBins).
  addBarrelEndcapSubsystem(builder, tracker, "LongStrips", "ls",
                           ActsPlugins::DD4hep::detail::kLongStripLayerFilter,
                           InnerCylinder, kLongStripBarrelZBins,
                           kLongStripEndcapRBins,
                           /*outerBoundaryZBins=*/kLongStripOuterBoundaryZBins);

  // Solenoid: a passive aluminum tube-shaped element outside LongStrips,
  // structurally identical to the beampipe/PST case. Unlike PST, its ODD XML
  // carries both a `layer_material surface="representing"` (-> outer face,
  // same convention as PST/beampipe) and a `boundary_material surface="inner"`
  // (-> inner face) on the same element.
  addPassiveCylinder(builder, tracker, "Solenoid",
                     {InnerCylinder, OuterCylinder}, kSolenoidZBins);

  addCaloMaterialCollector(builder, tracker);

  return finalizeOuter(builder, root, std::move(outerNode),
                       std::move(trackerNode), gctx, logger);
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

  auto outerNode = std::make_shared<CylinderContainerBlueprintNode>(
      "OpenDataDetector", AxisR);
  outerNode->setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
  auto& outer = *outerNode;

  addBeampipe(builder, outer);

  auto trackerNode = makeTrackerNode();
  auto& tracker = *trackerNode;

  using enum Face;
  addDirectLayerSubsystem(builder, tracker, "Pixels", "pix",
                          ActsPlugins::DD4hep::detail::kPixelLayerFilter,
                          OuterCylinder, kPixelBarrelZBins, kPixelEndcapRBins,
                          /*outerBoundaryZBins=*/std::nullopt);
  addPassiveCylinder(builder, tracker, "PST", {OuterCylinder}, kPstZBins);
  addDirectLayerSubsystem(builder, tracker, "ShortStrips", "ss",
                          ActsPlugins::DD4hep::detail::kShortStripLayerFilter,
                          InnerCylinder, kShortStripBarrelZBins,
                          kShortStripEndcapRBins,
                          /*outerBoundaryZBins=*/kShortStripBarrelZBins);
  addDirectLayerSubsystem(builder, tracker, "LongStrips", "ls",
                          ActsPlugins::DD4hep::detail::kLongStripLayerFilter,
                          InnerCylinder, kLongStripBarrelZBins,
                          kLongStripEndcapRBins,
                          /*outerBoundaryZBins=*/kLongStripOuterBoundaryZBins);
  addPassiveCylinder(builder, tracker, "Solenoid",
                     {InnerCylinder, OuterCylinder}, kSolenoidZBins);

  addCaloMaterialCollector(builder, tracker);

  return finalizeOuter(builder, root, std::move(outerNode),
                       std::move(trackerNode), gctx, logger);
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

  auto outerNode = std::make_shared<CylinderContainerBlueprintNode>(
      "OpenDataDetector", AxisR);
  outerNode->setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
  auto& outer = *outerNode;

  addBeampipe(builder, outer);

  auto trackerNode = makeTrackerNode();
  auto& tracker = *trackerNode;

  using enum Face;
  addDirectLayerGroupedSubsystem(builder, tracker, "Pixels", "pix",
                                 ActsPlugins::DD4hep::detail::kPixelLayerFilter,
                                 OuterCylinder, kPixelBarrelZBins,
                                 kPixelEndcapRBins,
                                 /*outerBoundaryZBins=*/std::nullopt);
  addPassiveCylinder(builder, tracker, "PST", {OuterCylinder}, kPstZBins);
  addDirectLayerGroupedSubsystem(
      builder, tracker, "ShortStrips", "ss",
      ActsPlugins::DD4hep::detail::kShortStripLayerFilter, InnerCylinder,
      kShortStripBarrelZBins, kShortStripEndcapRBins,
      /*outerBoundaryZBins=*/kShortStripBarrelZBins);
  addDirectLayerGroupedSubsystem(
      builder, tracker, "LongStrips", "ls",
      ActsPlugins::DD4hep::detail::kLongStripLayerFilter, InnerCylinder,
      kLongStripBarrelZBins, kLongStripEndcapRBins,
      /*outerBoundaryZBins=*/kLongStripOuterBoundaryZBins);
  addPassiveCylinder(builder, tracker, "Solenoid",
                     {InnerCylinder, OuterCylinder}, kSolenoidZBins);

  addCaloMaterialCollector(builder, tracker);

  return finalizeOuter(builder, root, std::move(outerNode),
                       std::move(trackerNode), gctx, logger);
}

}  // namespace ActsPlugins::DD4hep
