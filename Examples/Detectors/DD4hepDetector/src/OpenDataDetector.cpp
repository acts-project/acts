// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/OpenDataDetector.hpp"

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

void OpenDataDetector::construct(const Acts::GeometryContext& gctx) {
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

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
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

  outer.addCylinderContainer("Pixel", AxisZ, [&](auto& pixel) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});
    auto barrel =
        builder.layerHelper()
            .barrel()
            .setAxes("XYZ")
            .setPattern("PixelLayer\\d")
            .setContainer("PixelBarrel")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& elem, auto& layer) {
              layer.setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()
                      .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                          .layerType = SrfArrayNavPol::LayerType::Cylinder,
                          .bins = makeBinningFromConstants(
                              elem, std::regex{"PixelLayer(\\d+)"},
                              "pix_b{}_sf_b_phi", "pix_b_sf_b_z")})
                      .asUniquePtr());
            })
            .build();
    barrel->setAttachmentStrategy(AttachmentStrategy::First);

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {constant("pix_e_sf_b_r"), constant("pix_e_sf_b_phi")}})
            .asUniquePtr();

    auto negEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("PixelEndcapN\\d")
            .setContainer("PixelEndcapN")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& /*elem*/, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    negEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    auto posEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("PixelEndcapP\\d")
            .setContainer("PixelEndcapP")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& /*elem*/, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    posEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    pixel.addChild(barrel);
    pixel.addChild(negEndcap);
    pixel.addChild(posEndcap);
  });

  outer.addCylinderContainer("ShortStrip", AxisZ, [&](auto& sstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    auto barrel =
        builder.layerHelper()
            .barrel()
            .setAxes("XYZ")
            .setPattern("ShortStripLayer\\d")
            .setContainer("ShortStripBarrel")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& elem, auto& layer) {
              layer.setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()
                      .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                          .layerType = SrfArrayNavPol::LayerType::Cylinder,
                          .bins = makeBinningFromConstants(
                              elem, std::regex{"ShortStripLayer(\\d)"},
                              "ss_b{}_sf_b_phi", "ss_b_sf_b_z")})
                      .asUniquePtr());
            })
            .build();
    barrel->setAttachmentStrategy(AttachmentStrategy::First);

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {constant("ss_e_sf_b_r"), constant("ss_e_sf_b_phi")}})
            .asUniquePtr();

    auto posEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("ShortStripEndcapP\\d")
            .setContainer("ShortStripEndcapP")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& /*elem*/, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    posEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    auto negEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("ShortStripEndcapN\\d")
            .setContainer("ShortStripEndcapN")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& /*elem*/, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    negEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    sstrip.addChild(barrel);
    sstrip.addChild(posEndcap);
    sstrip.addChild(negEndcap);
  });

  outer.addCylinderContainer("LongStrip", AxisZ, [&](auto& lstrip) {
    auto envelope =
        ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm});

    auto barrel =
        builder.layerHelper()
            .barrel()
            .setAxes("XYZ")
            .setPattern("LongStripLayer\\d")
            .setContainer("LongStripBarrel")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& elem, auto& layer) {
              layer.setNavigationPolicyFactory(
                  NavigationPolicyFactory{}
                      .add<CylinderNavigationPolicy>()

                      .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                          .layerType = SrfArrayNavPol::LayerType::Cylinder,
                          .bins = makeBinningFromConstants(
                              elem, std::regex{"LongStripLayer(\\d)"},
                              "ls_b{}_sf_b_phi", "ls_b_sf_b_z")})
                      .asUniquePtr());
            })
            .build();
    barrel->setAttachmentStrategy(AttachmentStrategy::First);

    std::shared_ptr endcapPolicyFactory =
        NavigationPolicyFactory{}
            .add<CylinderNavigationPolicy>()
            .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
                .layerType = SrfArrayNavPol::LayerType::Disc,
                .bins = {constant("ls_e_sf_b_r"), constant("ls_e_sf_b_phi")}})
            .asUniquePtr();

    auto posEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("LongStripEndcapP\\d")
            .setContainer("LongStripEndcapP")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& /*elem*/, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    posEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    auto negEndcap =
        builder.layerHelper()
            .endcap()
            .setAxes("XZY")
            .setPattern("LongStripEndcapN\\d")
            .setContainer("LongStripEndcapN")
            .setEnvelope(envelope)
            .customize([&](const dd4hep::DetElement& /*elem*/, auto& layer) {
              layer.setNavigationPolicyFactory(endcapPolicyFactory);
            })
            .build();
    negEndcap->setAttachmentStrategy(AttachmentStrategy::First);

    lstrip.addChild(barrel);
    lstrip.addChild(posEndcap);
    lstrip.addChild(negEndcap);
  });

  // @TODO: Add plugin way to take this from xml

  BlueprintOptions options;

  m_trackingGeometry = root.construct(options, gctx, logger());
}

}  // namespace ActsExamples
