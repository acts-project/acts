// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/ePICDetector.hpp"

#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include <Acts/Geometry/Blueprint.hpp>
#include <Acts/Geometry/BlueprintOptions.hpp>
#include <Acts/Geometry/BoundarySurfaceFace.hpp>
#include <Acts/Geometry/ContainerBlueprintNode.hpp>
#include <Acts/Geometry/CylinderVolumeBounds.hpp>
#include <Acts/Geometry/Extent.hpp>
#include <Acts/Geometry/Layer.hpp>
#include <Acts/Navigation/CylinderNavigationPolicy.hpp>
#include <Acts/Navigation/SurfaceArrayNavigationPolicy.hpp>
#include <Acts/Navigation/TryAllNavigationPolicy.hpp>
#include <Acts/Surfaces/SurfaceArray.hpp>
#include <Acts/Utilities/AxisDefinitions.hpp>
#include <ActsPlugins/DD4hep/DD4hepConversionHelpers.hpp>

#include <algorithm>
#include <iterator>
#include <utility>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>

namespace ActsExamples {

ePICDetector::ePICDetector(const Config& cfg,
                           const Acts::GeometryContext& gctx)
    : DD4hepDetectorBase{cfg}, m_cfg{cfg} {
  ACTS_INFO("ePICDetector construct");
  construct(gctx);
}

auto ePICDetector::config() const -> const Config& {
  return m_cfg;
}

std::shared_ptr<ActsPlugins::DD4hepDetectorElement>
ePICDetector::defaultDetectorElementFactory(
    const dd4hep::DetElement& element, const std::string& axes, double scale) {
  return std::make_shared<ActsPlugins::DD4hepDetectorElement>(element, axes,
                                                              scale);
}

namespace {

std::optional<int> parseLayerNumber(const dd4hep::DetElement& elem,
                                    const std::regex& pattern) {
  std::cmatch match;

  if (!std::regex_match(elem.name(), match, pattern)) {
    return std::nullopt;
  }

  int n = std::stoi(match[1]);
  return n;
}

}  // namespace

void ePICDetector::construct(const Acts::GeometryContext& gctx) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using namespace Acts::UnitLiterals;
  using enum AxisDirection;

  ActsPlugins::DD4hep::BlueprintBuilder builder{
      {
          .dd4hepDetector = &dd4hepDetector(),
          .lengthScale = Acts::UnitConstants::cm,
      },
      logger().cloneWithSuffix("BlpBld")};

  // BARREL: XYZ
  // ENDCAP: XZY

  Blueprint::Config cfg;
  cfg.envelope[AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisR] = {0_mm, 20_mm};
  Blueprint root{cfg};

  // TODO DetElement acts_beampipe_central is beampipe
  // constants:
  // - IPBeampipeID,
  // - IPBeampipeUpstreamStraightLength,
  // - IPBeampipeDownstreamStraightLength
  auto& outer = root.addCylinderContainer("ePICDetector", AxisR);
  outer.addStaticVolume(
      Transform3::Identity(),
      std::make_unique<CylinderVolumeBounds>(0_mm, 20_mm, 1000_mm), "Beampipe");

  using AttachmentStrategy = Acts::VolumeAttachmentStrategy;
  using SrfArrayNavPol = Acts::SurfaceArrayNavigationPolicy;

  auto constant = [this](const std::string& name) -> int {
    return dd4hepDetector().constant<int>(name);
  };

  auto makeBinningFromConstants =
      [&](const dd4hep::DetElement& elem, const std::regex& pattern,
          const std::string& constant0, const std::string& constant1) {
        std::cmatch match;

        if (!std::regex_match(elem.name(), match, pattern)) {
          throw std::runtime_error(std::format(
              "Could not extract layer number from {}", elem.name()));
        }

        if (auto n = parseLayerNumber(elem, pattern)) {
          return std::pair{
              constant(std::vformat(constant0, std::make_format_args(*n))),
              constant(std::vformat(constant1, std::make_format_args(*n)))};
        } else {
          throw std::runtime_error(std::format(
              "Could not extract layer number from {}", elem.name()));
        }
      };

  // VertexBarrel
  // TODO common pattern: nphi by variable, nz by fixed constant
  // which may need makeBinningFromConstants with variant<string,int>
  outer.addCylinderContainer("VertexBarrel", AxisZ, [&](auto& pixel) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});
    auto barrel =
        builder.layerHelper()
            .barrel()
            .setAxes("XYZ")
            .setPattern("VertexBarrel_layer\\d")
            .setContainer("VertexBarrel")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& elem, auto& layer) {
              layer.setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()
                      .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                          .layerType = SrfArrayNavPol::LayerType::Cylinder,
                          .bins = makeBinningFromConstants(
                              elem, std::regex{"VertexBarrel_layer(\\d)"},
                              "VertexBarrelL{}_nphi", "VertexBarrelL{}_nz")})
                      .asUniquePtr());
            })
            .build();
    barrel->setAttachmentStrategy(AttachmentStrategy::First);

    pixel.addChild(barrel);
  });

  // SagittaSiBarrel
  outer.addCylinderContainer("SagittaSiBarrel", AxisZ, [&](auto& sstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    auto barrel =
        builder.layerHelper()
            .barrel()
            .setAxes("XYZ")
            .setPattern("SagittaSiBarrel_layer\\d")
            .setContainer("SagittaSiBarrel")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& elem, auto& layer) {
              layer.setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()
                      .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                          .layerType = SrfArrayNavPol::LayerType::Cylinder,
                          .bins = {constant("SiBarrelStave1_count"), 100}})
                      .asUniquePtr());
            })
            .build();
    barrel->setAttachmentStrategy(AttachmentStrategy::First);

    sstrip.addChild(barrel);
  });

  // OuterSiBarrel
  outer.addCylinderContainer("OuterSiBarrel", AxisZ, [&](auto& sstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    auto barrel =
        builder.layerHelper()
            .barrel()
            .setAxes("XYZ")
            .setPattern("OuterSiBarrel_layer\\d")
            .setContainer("OuterSiBarrel")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& elem, auto& layer) {
              layer.setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()
                      .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                          .layerType = SrfArrayNavPol::LayerType::Cylinder,
                          .bins = {128, 100}})
                      .asUniquePtr());
            })
            .build();
    barrel->setAttachmentStrategy(AttachmentStrategy::First);

    sstrip.addChild(barrel);
  });

  // InnerTrackerEndcapP
  outer.addCylinderContainer("InnerTrackerEndcapP", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {5*constant("SiTrackerEndcapMod_count"), 100}})
            .asUniquePtr();

    auto posEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("InnerTrackerEndcapP_layer\\d_P")
            .setContainer("InnerTrackerEndcapP")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement&, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    posEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(posEndcap);
  });

  // InnerTrackerEndcapN
  outer.addCylinderContainer("InnerTrackerEndcapN", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});
        
    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {5*constant("SiTrackerEndcapMod_count"), 100}})
            .asUniquePtr();
        
    auto negEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("InnerTrackerEndcapN_layer\\d_N")
            .setContainer("InnerTrackerEndcapN")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement&, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    negEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(negEndcap);
  });

  // MiddleTrackerEndcapP
  outer.addCylinderContainer("MiddleTrackerEndcapP", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {5*constant("SiTrackerEndcapMod_count"), 100}})
            .asUniquePtr();

    auto posEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("MiddleTrackerEndcapP_layer\\d_P")
            .setContainer("MiddleTrackerEndcapP")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement&, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    posEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(posEndcap);
  });

  // MiddleTrackerEndcapN
  outer.addCylinderContainer("MiddleTrackerEndcapN", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {5*constant("SiTrackerEndcapMod_count"), 100}})
            .asUniquePtr();

    auto negEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("MiddleTrackerEndcapN_layer\\d_N")
            .setContainer("MiddleTrackerEndcapN")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement&, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    negEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(negEndcap);
  });

  // OuterTrackerEndcapP
  outer.addCylinderContainer("OuterTrackerEndcapP", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {5*constant("SiTrackerEndcapMod_count"), 100}})
            .asUniquePtr();

    auto posEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("OuterTrackerEndcapP_layer\\d_P")
            .setContainer("OuterTrackerEndcapP")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement&, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    posEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(posEndcap);
  });

  // OuterTrackerEndcapN
  outer.addCylinderContainer("OuterTrackerEndcapN", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {5*constant("SiTrackerEndcapMod_count"), 100}})
            .asUniquePtr();

    auto negEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("OuterTrackerEndcapN_layer\\d_N")
            .setContainer("OuterTrackerEndcapN")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement&, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    negEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(negEndcap);
  });


  // InnerMPGDBarrel
  outer.addCylinderContainer("InnerMPGDBarrel", AxisZ, [&](auto& sstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    auto barrel =
        builder.layerHelper()
            .barrel()
            .setAxes("XYZ")
            .setPattern("InnerMPGDBarrel_layer\\d")
            .setContainer("InnerMPGDBarrel")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& elem, auto& layer) {
              layer.setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()
                      .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                          .layerType = SrfArrayNavPol::LayerType::Cylinder,
                          .bins = {24*constant("MPGDBarrelStave_count"), 100}})
                      .asUniquePtr());
            })                                                
            .build();
    barrel->setAttachmentStrategy(AttachmentStrategy::First);

    sstrip.addChild(barrel);
  });

  // MPGDOuterBarrel
  outer.addCylinderContainer("MPGDOuterBarrel", AxisZ, [&](auto& sstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    auto barrel =
        builder.layerHelper()
            .barrel()
            .setAxes("XYZ")
            .setPattern("MPGDOuterBarrel_layer\\d")
            .setContainer("MPGDOuterBarrel")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& elem, auto& layer) {
              layer.setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()
                      .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                          .layerType = SrfArrayNavPol::LayerType::Cylinder,
                          .bins = {10*constant("MPGDOuterBarrelModule_count"), 100}})
                      .asUniquePtr());
            })
            .build();
    barrel->setAttachmentStrategy(AttachmentStrategy::First);

    sstrip.addChild(barrel);
  });

  // ForwardMPGD
  outer.addCylinderContainer("ForwardMPGD", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {constant("ForwardMPGDEndcapMod_count"), 30}})
            .asUniquePtr();

    auto posEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("ForwardMPGD_layer\\d_P")
            .setContainer("ForwardMPGD")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement&, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    posEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(posEndcap);
  });

  // BackwardMPGD
  outer.addCylinderContainer("BackwardMPGD", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {constant("BackwardMPGDEndcapMod_count"), 50}})
            .asUniquePtr();

    auto negEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("BackwardMPGD_layer\\d_N")
            .setContainer("BackwardMPGD")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement&, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    negEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(negEndcap);
  });

  // BarrelTOF
  outer.addCylinderContainer("BarrelTOF", AxisZ, [&](auto& sstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    auto barrel =
        builder.layerHelper()
            .barrel()
            .setAxes("XYZ")
            .setPattern("BarrelTOF_layer\\d")
            .setContainer("BarrelTOF")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& elem, auto& layer) {
              layer.setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()
                      .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                          .layerType = SrfArrayNavPol::LayerType::Cylinder,
                          .bins = {constant("BarrelTOF_Module_nphi"), 100}})
                      .asUniquePtr());
            })
            .build();
    barrel->setAttachmentStrategy(AttachmentStrategy::First); 
    
    sstrip.addChild(barrel);
  });

  // ForwardTOF
  outer.addCylinderContainer("ForwardTOF", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {30, 30}})
            .asUniquePtr();

    auto posEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("ForwardTOF_layer1")
            //.setPattern("ForwardTOF_layer\\d")
            // ^ FIXME LayerBlueprintNode: no surfaces provided for ForwardTOF_layer2
            .setContainer("ForwardTOF")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement&, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    posEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(posEndcap);
  });

  // B0 Tracker
  // FIXME VolumeStack requires at least one volume
  /*
  outer.addCylinderContainer("B0Tracker", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {2*constant("B0TrackerLayerSmallMod_nModules"), 12}})
            .asUniquePtr();

    auto posEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("B0Tracker_layer\\d")
            .setContainer("B0Tracker")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement&, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    posEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(posEndcap);     
  });
  */

  ACTS_INFO("ePICDetector blueprint");

  // @TODO: Add plugin way to take this from xml

  BlueprintOptions options;

  m_trackingGeometry = root.construct(options, gctx, logger());
}

}  // namespace ActsExamples
