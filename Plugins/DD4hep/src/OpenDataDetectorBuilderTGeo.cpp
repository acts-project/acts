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
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "ActsPlugins/DD4hep/OpenDataDetectorBuilder.hpp"
#include "ActsPlugins/Root/BlueprintBuilder.hpp"
#include "ActsPlugins/Root/TGeoSurfaceConverter.hpp"

#include <array>
#include <format>
#include <memory>
#include <optional>
#include <regex>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>

#include <DD4hep/Detector.h>

#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"

namespace ActsPlugins::DD4hep {

namespace {

const std::regex kTGeoPixelLayerFilter{"(?:PixelBarrel|PixelEndcap[NP])(\\d)"};
const std::regex kTGeoShortStripLayerFilter{
    "(?:ShortStripBarrel|ShortStripEndcap[NP])(\\d)"};
const std::regex kTGeoLongStripLayerFilter{
    "(?:LongStripBarrel|LongStripEndcap[NP])(\\d)"};

const std::regex kTGeoPixelBarrelLayerFilter{"PixelBarrel\\d"};
const std::regex kTGeoShortStripBarrelLayerFilter{"ShortStripBarrel\\d"};
const std::regex kTGeoLongStripBarrelLayerFilter{"LongStripBarrel\\d"};
constexpr std::string_view kTGeoBeampipeName = "BeamPipe";

struct TGeoLayerBinning {
  std::span<const std::size_t> barrelPhiBins = {};
  std::size_t barrelZBins = 0u;
  std::array<std::size_t, 2> endcapBins = {0u, 0u};

  std::array<std::size_t, 2> binsFor(bool isBarrelLayer, int layerIdx) const {
    if (!isBarrelLayer) {
      return endcapBins;
    }

    if (layerIdx < 0 ||
        static_cast<std::size_t>(layerIdx) >= barrelPhiBins.size()) {
      throw std::runtime_error(std::format(
          "No TGeo barrel binning configured for layer {}", layerIdx));
    }

    return {barrelPhiBins[static_cast<std::size_t>(layerIdx)], barrelZBins};
  }
};

constexpr std::array<std::size_t, 4> kPixelBarrelPhiBins = {16u, 32u, 52u, 78u};
constexpr TGeoLayerBinning kPixelLayerBinning{
    .barrelPhiBins = kPixelBarrelPhiBins,
    .barrelZBins = 14u,
    .endcapBins = {2u, 36u},
};

constexpr std::array<std::size_t, 4> kShortStripBarrelPhiBins = {40u, 56u, 78u,
                                                                 102u};
constexpr TGeoLayerBinning kShortStripLayerBinning{
    .barrelPhiBins = kShortStripBarrelPhiBins,
    .barrelZBins = 21u,
    .endcapBins = {3u, 42u},
};

constexpr std::array<std::size_t, 2> kLongStripBarrelPhiBins = {60u, 80u};
constexpr TGeoLayerBinning kLongStripLayerBinning{
    .barrelPhiBins = kLongStripBarrelPhiBins,
    .barrelZBins = 21u,
    .endcapBins = {2u, 48u},
};

bool hasSensitiveMaterial(const TGeoBlueprintBuilderBackend::Element& element) {
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

auto makeTGeoLayerCustomizer(const BlueprintBuilder& builder,
                             const TGeoLayerBinning& binning,
                             const std::regex& layerFilter) {
  return [&builder, &binning, layerFilter](
             const std::optional<TGeoBlueprintBuilderBackend::Element>& elem,
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
    navCfg.layerType = isBarrelLayer ? Cylinder : Disc;

    const auto bins = binning.binsFor(isBarrelLayer, layerIdx);
    navCfg.bins = {bins[0], bins[1]};

    layer.setNavigationPolicyFactory(Acts::NavigationPolicyFactory{}
                                         .add<Acts::CylinderNavigationPolicy>()
                                         .add<SrfArrayNavPol>(navCfg)
                                         .asUniquePtr());
  };
}

std::shared_ptr<Acts::Experimental::StaticBlueprintNode> makeTGeoBeampipeNode(
    const BlueprintBuilder& builder) {
  const auto beampipeElement =
      builder.findDetElementByName(std::string{kTGeoBeampipeName});
  if (!beampipeElement.has_value()) {
    throw std::runtime_error(
        std::format("Could not find beampipe element '{}'", kTGeoBeampipeName));
  }

  const auto& node = builder.backend().nodeOf(*beampipeElement);
  const auto tgTransform = builder.backend().transformOf(*beampipeElement);
  auto [bounds, transform, thickness] =
      ActsPlugins::TGeoSurfaceConverter::cylinderComponents(
          *node.GetVolume()->GetShape(), tgTransform.GetRotationMatrix(),
          tgTransform.GetTranslation(), "XYZ", Acts::UnitConstants::cm);
  (void)thickness;

  if (bounds == nullptr) {
    throw std::runtime_error(
        "Beampipe element shape could not be converted to cylinder.");
  }

  auto volumeBounds = std::make_shared<Acts::CylinderVolumeBounds>(
      0., bounds->get(Acts::CylinderBounds::eR),
      bounds->get(Acts::CylinderBounds::eHalfLengthZ));
  auto volume = std::make_unique<Acts::TrackingVolume>(
      transform, volumeBounds, std::string{kTGeoBeampipeName});
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
  std::string_view barrelName;
  std::string_view negativeEndcapName;
  std::string_view positiveEndcapName;
  const TGeoLayerBinning& binning;
  const std::regex& layerFilter;
  const std::regex& barrelLayerFilter;
  const std::regex& negativeEndcapLayerFilter;
  const std::regex& positiveEndcapLayerFilter;
};

void addTGeoSubsystem(const BlueprintBuilder& builder,
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
                        .onLayer(makeTGeoLayerCustomizer(builder, spec.binning,
                                                         spec.layerFilter))
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
          .onLayer(
              makeTGeoLayerCustomizer(builder, spec.binning, spec.layerFilter))
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
          .onLayer(
              makeTGeoLayerCustomizer(builder, spec.binning, spec.layerFilter))
          .build();
  configureSubsystemNode(*positiveEndcapNode);
  subsystemNode->addChild(std::move(positiveEndcapNode));

  parent.addChild(std::move(subsystemNode));
}

}  // namespace

std::unique_ptr<Acts::TrackingGeometry>
buildOpenDataDetectorBarrelEndcapViaTGeo(const TGeoNode& rootNode,
                                         const Acts::GeometryContext& gctx,
                                         const Acts::Logger& logger) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using enum AxisDirection;

  TGeoBlueprintBuilderBackend::Config cfg;
  cfg.root = &rootNode;
  cfg.lengthScale = Acts::UnitConstants::cm;
  cfg.sensitivePredicate = hasSensitiveMaterial;

  BlueprintBuilder builder{cfg, logger.cloneWithSuffix("TGeoBlpBld")};

  Blueprint::Config blueprintCfg;
  blueprintCfg.envelope = ActsPlugins::DD4hep::detail::kBlueprintEnvelope;
  Blueprint root{blueprintCfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  outer.addChild(makeTGeoBeampipeNode(builder));
  addTGeoSubsystem(
      builder, outer,
      {.assembly = "Pixels",
       .barrelName = "PixelBarrel",
       .negativeEndcapName = "PixelEndcapN",
       .positiveEndcapName = "PixelEndcapP",
       .binning = kPixelLayerBinning,
       .layerFilter = kTGeoPixelLayerFilter,
       .barrelLayerFilter = kTGeoPixelBarrelLayerFilter,
       .negativeEndcapLayerFilter = detail::kPixelNegativeEndcapLayerFilter,
       .positiveEndcapLayerFilter = detail::kPixelPositiveEndcapLayerFilter});
  addTGeoSubsystem(builder, outer,
                   {.assembly = "ShortStrips",
                    .barrelName = "ShortStripBarrel",
                    .negativeEndcapName = "ShortStripEndcapN",
                    .positiveEndcapName = "ShortStripEndcapP",
                    .binning = kShortStripLayerBinning,
                    .layerFilter = kTGeoShortStripLayerFilter,
                    .barrelLayerFilter = kTGeoShortStripBarrelLayerFilter,
                    .negativeEndcapLayerFilter =
                        detail::kShortStripNegativeEndcapLayerFilter,
                    .positiveEndcapLayerFilter =
                        detail::kShortStripPositiveEndcapLayerFilter});
  addTGeoSubsystem(
      builder, outer,
      {.assembly = "LongStrips",
       .barrelName = "LongStripBarrel",
       .negativeEndcapName = "LongStripEndcapN",
       .positiveEndcapName = "LongStripEndcapP",
       .binning = kLongStripLayerBinning,
       .layerFilter = kTGeoLongStripLayerFilter,
       .barrelLayerFilter = kTGeoLongStripBarrelLayerFilter,
       .negativeEndcapLayerFilter = detail::kLongStripNegativeEndcapLayerFilter,
       .positiveEndcapLayerFilter =
           detail::kLongStripPositiveEndcapLayerFilter});

  return root.construct(BlueprintOptions{}, gctx, logger);
}

}  // namespace ActsPlugins::DD4hep
