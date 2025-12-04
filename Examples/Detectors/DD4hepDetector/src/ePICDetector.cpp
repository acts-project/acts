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

  //
  // DEFINE DETECTORS
  //

  // VertexBarrel
  auto VertexBarrel =
      builder.layerHelper()
          .barrel()
          .setAxes("XYZ")
          .setPattern("VertexBarrel_layer\\d")
          .setContainer("VertexBarrel")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
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
  VertexBarrel->setAttachmentStrategy(AttachmentStrategy::First);

  // SagittaSiBarrel
  // FIXME Volumes are not aligned: translation in x or y
  // Requires changing number of modules to multiple of 4
  auto SagittaSiBarrel =
      builder.layerHelper()
          .barrel()
          .setAxes("XYZ")
          .setPattern("SagittaSiBarrel_layer\\d")
          .setContainer("SagittaSiBarrel")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
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
  SagittaSiBarrel->setAttachmentStrategy(AttachmentStrategy::First);

  // OuterSiBarrel
  // FIXME Volumes are not aligned: translation in x or y
  // Requires changing number of modules to multiple of 4
  auto OuterSiBarrel =
      builder.layerHelper()
          .barrel()
          .setAxes("XYZ")
          .setPattern("OuterSiBarrel_layer\\d")
          .setContainer("OuterSiBarrel")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
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
  OuterSiBarrel->setAttachmentStrategy(AttachmentStrategy::First);

  // endcapPolicyFactory
  std::shared_ptr SiTrackerEndcapPolicyFactory =
    NavigationPolicyFactory{}
        .add<CylinderNavigationPolicy>()
        .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
            .layerType = SrfArrayNavPol::LayerType::Disc,
            .bins = {5*constant("SiTrackerEndcapMod_count"), 100}})
        .asUniquePtr();

  // InnerTrackerEndcapP
  auto InnerTrackerEndcapP =
      builder.layerHelper()
          .endcap()
          .setAxes("XZY")
          .setPattern("InnerTrackerEndcapP_layer\\d_P")
          .setContainer("InnerTrackerEndcapP")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
          .customize([&](const dd4hep::DetElement&, auto& layer) {
            layer.setNavigationPolicyFactory(SiTrackerEndcapPolicyFactory);
          })
          .build();
  InnerTrackerEndcapP->setAttachmentStrategy(AttachmentStrategy::First);

  // InnerTrackerEndcapN
  auto InnerTrackerEndcapN =
      builder.layerHelper()
          .endcap()
          .setAxes("XZY")
          .setPattern("InnerTrackerEndcapN_layer\\d_N")
          .setContainer("InnerTrackerEndcapN")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
          .customize([&](const dd4hep::DetElement&, auto& layer) {
            layer.setNavigationPolicyFactory(SiTrackerEndcapPolicyFactory);
          })
          .build();
  InnerTrackerEndcapN->setAttachmentStrategy(AttachmentStrategy::First);

  // MiddleTrackerEndcapP
  auto MiddleTrackerEndcapP =
      builder.layerHelper()
          .endcap()
          .setAxes("XZY")
          .setPattern("MiddleTrackerEndcapP_layer\\d_P")
          .setContainer("MiddleTrackerEndcapP")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
          .customize([&](const dd4hep::DetElement&, auto& layer) {
            layer.setNavigationPolicyFactory(SiTrackerEndcapPolicyFactory);
          })
          .build();
  MiddleTrackerEndcapP->setAttachmentStrategy(AttachmentStrategy::First);

  // MiddleTrackerEndcapN
  auto MiddleTrackerEndcapN =
      builder.layerHelper()
          .endcap()
          .setAxes("XZY")
          .setPattern("MiddleTrackerEndcapN_layer\\d_N")
          .setContainer("MiddleTrackerEndcapN")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
          .customize([&](const dd4hep::DetElement&, auto& layer) {
            layer.setNavigationPolicyFactory(SiTrackerEndcapPolicyFactory);
          })
          .build();
  MiddleTrackerEndcapN->setAttachmentStrategy(AttachmentStrategy::First);

  // OuterTrackerEndcapP
  auto OuterTrackerEndcapP =
      builder.layerHelper()
          .endcap()
          .setAxes("XZY")
          .setPattern("OuterTrackerEndcapP_layer\\d_P")
          .setContainer("OuterTrackerEndcapP")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
          .customize([&](const dd4hep::DetElement&, auto& layer) {
            layer.setNavigationPolicyFactory(SiTrackerEndcapPolicyFactory);
          })
          .build();
  OuterTrackerEndcapP->setAttachmentStrategy(AttachmentStrategy::First);

  // OuterTrackerEndcapN
  auto OuterTrackerEndcapN =
      builder.layerHelper()
          .endcap()
          .setAxes("XZY")
          .setPattern("OuterTrackerEndcapN_layer\\d_N")
          .setContainer("OuterTrackerEndcapN")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
          .customize([&](const dd4hep::DetElement&, auto& layer) {
            layer.setNavigationPolicyFactory(SiTrackerEndcapPolicyFactory);
          })
          .build();
  OuterTrackerEndcapN->setAttachmentStrategy(AttachmentStrategy::First);

  // ForwardMPGD
  std::shared_ptr ForwardMPGDPolicyFactory =
      NavigationPolicyFactory{}
          .add<CylinderNavigationPolicy>()
          .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
              .layerType = SrfArrayNavPol::LayerType::Disc,
              .bins = {constant("ForwardMPGDEndcapMod_count"), 30}})
          .asUniquePtr();
  auto ForwardMPGD =
      builder.layerHelper()
          .endcap()
          .setAxes("XZY")
          .setPattern("ForwardMPGD_layer\\d_P")
          .setContainer("ForwardMPGD")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
          .customize([&](const dd4hep::DetElement&, auto& layer) {
            layer.setNavigationPolicyFactory(ForwardMPGDPolicyFactory);
          })
          .build();
  ForwardMPGD->setAttachmentStrategy(AttachmentStrategy::First);

  // BackwardMPGD
  std::shared_ptr BackwardMPGDPolicyFactory =
      NavigationPolicyFactory{}
          .add<CylinderNavigationPolicy>()
          .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
              .layerType = SrfArrayNavPol::LayerType::Disc,
              .bins = {constant("BackwardMPGDEndcapMod_count"), 50}})
          .asUniquePtr();
  auto BackwardMPGD =
      builder.layerHelper()
          .endcap()
          .setAxes("XZY")
          .setPattern("BackwardMPGD_layer\\d_N")
          .setContainer("BackwardMPGD")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
          .customize([&](const dd4hep::DetElement&, auto& layer) {
            layer.setNavigationPolicyFactory(BackwardMPGDPolicyFactory);
          })
          .build();
  BackwardMPGD->setAttachmentStrategy(AttachmentStrategy::First);

  // InnerMPGDBarrel
  auto InnerMPGDBarrel =
      builder.layerHelper()
          .barrel()
          .setAxes("XYZ")
          .setPattern("InnerMPGDBarrel_layer\\d")
          .setContainer("InnerMPGDBarrel")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
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
  InnerMPGDBarrel->setAttachmentStrategy(AttachmentStrategy::First);

  // BarrelTOF
  auto BarrelTOF =
      builder.layerHelper()
          .barrel()
          .setAxes("XYZ")
          .setPattern("BarrelTOF_layer\\d")
          .setContainer("BarrelTOF")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
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
  BarrelTOF->setAttachmentStrategy(AttachmentStrategy::First);

  // MPGDOuterBarrel
  auto MPGDOuterBarrel =
      builder.layerHelper()
          .barrel()
          .setAxes("XYZ")
          .setPattern("MPGDOuterBarrel_layer\\d")
          .setContainer("MPGDOuterBarrel")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
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
  MPGDOuterBarrel->setAttachmentStrategy(AttachmentStrategy::First);

  // ForwardTOF
  // FIXME Volumes are not aligned: translation in x or y
  std::shared_ptr ForwardTOFPolicyFactory =
      NavigationPolicyFactory{}
          .add<CylinderNavigationPolicy>()
          .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
              .layerType = SrfArrayNavPol::LayerType::Disc,
              .bins = {30, 30}})
          .asUniquePtr();
  auto ForwardTOF =
      builder.layerHelper()
          .endcap()
          .setAxes("XZY")
          .setPattern("ForwardTOF_layer1")
          //.setPattern("ForwardTOF_layer\\d")
          // ^ FIXME LayerBlueprintNode: no surfaces provided for ForwardTOF_layer2
          .setContainer("ForwardTOF")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}))
          .customize([&](const dd4hep::DetElement&, auto& layer) {
            layer.setNavigationPolicyFactory(ForwardTOFPolicyFactory);
          })
          .build();
  ForwardTOF->setAttachmentStrategy(AttachmentStrategy::First);

  // B0Tracker (OFF AXIS)
  // FIXME VolumeStack requires at least one volume
  /*
  std::shared_ptr B0TrackerPolicyFactory =
      NavigationPolicyFactory{}
          .add<CylinderNavigationPolicy>()
          .add<SrfArrayNavPol>(SrfArrayNavPol::Config{
              .layerType = SrfArrayNavPol::LayerType::Disc,
              .bins = {2*constant("B0TrackerLayerSmallMod_nModules"), 12}})
          .asUniquePtr();
  auto B0Tracker =
      builder.layerHelper()
          .endcap()
          .setAxes("XZY")
          .setPattern("B0Tracker_layer\\d")
          .setContainer("B0Tracker")
          .setEnvelope(
            ExtentEnvelope{}.set(AxisZ, {5_mm, 5_mm}).set(AxisR, {5_mm, 5_mm}
          .customize([&](const dd4hep::DetElement&, auto& layer) {
            layer.setNavigationPolicyFactory(B0TrackerPolicyFactory);
          })
          .build();
  B0Tracker->setAttachmentStrategy(AttachmentStrategy::First);
  */

  //
  // PLACE IN NESTED CONTAINERS
  //

  // TODO DetElement acts_beampipe_central is beampipe
  // constants:
  // - IPBeampipeID,
  // - IPBeampipeUpstreamStraightLength,
  // - IPBeampipeDownstreamStraightLength

  // Note: easiest to think from inside to outside
  root.addCylinderContainer("Tracker4", AxisZ, [&](auto& tracker4) {
    tracker4.addCylinderContainer("Tracker3", AxisR, [&](auto& tracker3) {
      tracker3.addStaticVolume(
          Transform3::Identity(),
          std::make_unique<CylinderVolumeBounds>(0_mm, 20_mm, 100_mm), "Beampipe");
      tracker3.addCylinderContainer("Tracker2", AxisZ, [&](auto& tracker2) {
        tracker2.addChild(BackwardMPGD);               // r=[65–405], z=[−1462,−1324]
        tracker2.addChild(OuterTrackerEndcapN);        // r=[32–426], z=[−1275,−895]
        tracker2.addChild(MiddleTrackerEndcapN);       // r=[32-420], z ~ -450
        tracker2.addCylinderContainer("Tracker1", AxisR, [&](auto& tracker1) {
          tracker1.addCylinderContainer("Tracker0", AxisZ, [&](auto& tracker0) {
            tracker0.addChild(InnerTrackerEndcapN);    // r=[32-245], z ~ -250
            tracker0.addChild(VertexBarrel);           // r=[33–130], z=[−135,+135]
            tracker0.addChild(InnerTrackerEndcapP);    // r=[32-245], z ~ +250
          });
          tracker1.addChild(SagittaSiBarrel);          // r=[258-275], z=[-256,+256]
          tracker1.addChild(OuterSiBarrel);            // r=[413–430], z=[−402,+402]
        });
        tracker2.addChild(MiddleTrackerEndcapP);       // r=[32-420], z ~ +450
        tracker2.addChild(OuterTrackerEndcapP);        // r=[34–426], z=[+695,+1355]
        tracker2.addChild(ForwardMPGD);                // r=[76–405], z=[+1249,+1387]
      });
      tracker3.addChild(InnerMPGDBarrel);              // r=[547–589], z=[−1192,+1192]
      tracker3.addChild(BarrelTOF);                    // r=[629–654], z=[−1285,+1285]
      tracker3.addChild(MPGDOuterBarrel);              // r=[731–762], z=[−1700,+1700]
    });
    tracker4.addChild(ForwardTOF);                     // r=[101–602], z ~ +1861
    //tracker4.addChild(B0Tracker);                      // r=[35-150], z=[+5895,+6705] OFF-AXIS
  });

  // @TODO: Add plugin way to take this from xml

  BlueprintOptions options;

  m_trackingGeometry = root.construct(options, gctx, logger());
}

}  // namespace ActsExamples
