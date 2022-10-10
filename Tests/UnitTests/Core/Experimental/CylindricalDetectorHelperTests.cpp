// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/CylindricalDetectorHelper.hpp"
#include "Acts/Experimental/Detector.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/detail/NavigationStateUpdators.hpp"
#include "Acts/Experimental/detail/PortalGenerators.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <exception>
#include <memory>

// Remove later - debug only
#include "Acts/Plugins/ActSVG/DetectorSvgConverter.hpp"

// Colorize in blue
actsvg::style::color red = actsvg::style::color({{255, 0, 0}});
actsvg::style::color green = actsvg::style::color({{0, 255, 0}});
actsvg::style::color blue = actsvg::style::color({{0, 0, 255}});
actsvg::style::color black = actsvg::style::color({{0, 0, 0}});
actsvg::style::color olive = actsvg::style::color({{128, 128, 0}});
actsvg::style::color orange = actsvg::style::color({{255, 153, 0}});
actsvg::style::color cyan = actsvg::style::color({{0, 255, 255}});
actsvg::style::color magenta = actsvg::style::color({{255, 0, 255}});
actsvg::style::color brown = actsvg::style::color({{153, 102, 51}});
actsvg::style::color purple = actsvg::style::color({{153, 0, 204}});
std::vector<actsvg::style::color> colors = {
    red, green, blue, black, olive, orange, cyan, magenta, brown, purple};

using namespace Acts::Experimental;

/// Trial and error volume finder
///
/// @param gctx is the geometry context of this call
/// @param detector is the detector
/// @param position is the position
///
/// @return a detector volume pointer or null
const DetectorVolume* trialAndError(const Acts::GeometryContext& gctx,
                                    const Detector& detector,
                                    const Acts::Vector3& position) {
  for (const auto v : detector.volumes()) {
    if (v->inside(gctx, position)) {
      return v;
    }
  }
  return nullptr;
}

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(ConnectPerfectlyFitInR) {
  // Test with different transforms
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};

  // Test with different opening angles
  std::vector<Acts::ActsScalar> openings = {M_PI};

  std::vector<Acts::ActsScalar> radii = {0., 10., 100., 200.};
  Acts::ActsScalar halfZ = 100;

  auto portalGenerator = detail::defaultPortalGenerator();
  auto navigationStateUpdator = detail::defaultPortalProvider();

  // This should work for generic transforms
  for (auto [it, tf] : Acts::enumerate(transforms)) {
    std::string trfStr = "transform_" + std::to_string(it);
    // This should work for full cylinder and sector openings
    for (auto [io, openings] : Acts::enumerate(openings)) {
      std::string opStr = "_opening_" + std::to_string(it);

      std::vector<std::shared_ptr<DetectorVolume>> rVolumes = {};
      // Create the voluems
      for (auto [i, r] : Acts::enumerate(radii)) {
        if (i > 0) {
          auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
              radii[i - 1u], r, halfZ);
          rVolumes.push_back(DetectorVolumeFactory::construct(
              portalGenerator, tContext, "Cylinder_r" + std::to_string(i), tf,
              std::move(cBounds), detail::defaultPortalProvider()));
        }
      }

      connectVolumesInR(tContext, rVolumes);

      // A detector construction that should work
      DetectorVolumeFinder trialAndErrorFinder;
      trialAndErrorFinder.connect<&trialAndError>();

      auto detector =
          Detector::makeShared("DetectorInR", rVolumes, trialAndErrorFinder);

      // Remove later
      // ------------------------------------------------------------
      Acts::Svg::Style portalStyle;
      portalStyle.strokeHighlights = {"mouseover", "mouseout"};
      portalStyle.strokeHighlightWidth = 3.;
      portalStyle.strokeHighlightColor = {255, 255, 0};

      Acts::Svg::DetectorConverter::Options detectorOptions;
      detectorOptions.volumeOptions.portalOptions.surfaceOptions.style =
          portalStyle;

      auto pDetector = Acts::Svg::DetectorConverter::convert(
          tContext, *detector, detectorOptions);
      pDetector._name = detector->name();

      for (auto& c : colors) {
        c._opacity = 0.1;
      }

      pDetector.colorize(colors);

      // As zr view
      auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
      Acts::Svg::toFile({dv_zr},
                        "DetectorVolumesInR_" + trfStr + opStr + "_zr.svg");

      auto dv_xy = Acts::Svg::View::xy(pDetector, pDetector._name);
      Acts::Svg::toFile({dv_xy},
                        "DetectorVolumesInR_" + trfStr + opStr + "_xy.svg");
    }
  }
}

BOOST_AUTO_TEST_CASE(PerfectFitAttachmentZ) {
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};
  std::vector<std::array<Acts::ActsScalar, 2>> radii = {{0., 100.},
                                                        {20., 120.}};
  std::vector<Acts::ActsScalar> zValues = {-100., -20, 10., 100., 200.};

  auto portalGenerator = detail::defaultPortalGenerator();
  auto navigationStateUpdator = detail::defaultPortalProvider();

  for (auto [it, t] : Acts::enumerate(transforms)) {
    std::string trfStr = "transform_" + std::to_string(it);
    for (auto [ir, r] : Acts::enumerate(radii)) {
      std::string radStr = "radii_" + std::to_string(ir);
      std::vector<std::shared_ptr<DetectorVolume>> zVolumes = {};
      for (auto [i, z] : Acts::enumerate(zValues)) {
        if (i > 0) {
          auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
              r[0], r[1], 0.5 * (z - zValues[i - 1u]));
          // z center
          Acts::ActsScalar zCenter = 0.5 * (z + zValues[i - 1u]);
          Acts::Transform3 ti = Acts::Transform3::Identity();
          ti.pretranslate(t.translation() +
                          zCenter * t.rotation().matrix().col(2));
          ti.prerotate(t.rotation());
          // create the volume
          zVolumes.push_back(DetectorVolumeFactory::construct(
              portalGenerator, tContext, "Cylinder_z" + std::to_string(i), ti,
              std::move(cBounds), detail::defaultPortalProvider()));
        }
      }
      // Now call the connector
      connectVolumesInZ(tContext, zVolumes);

      // A detector construction that should work
      DetectorVolumeFinder trialAndErrorFinder;
      trialAndErrorFinder.connect<&trialAndError>();

      auto detector =
          Detector::makeShared("DetectorInZ", zVolumes, trialAndErrorFinder);

      // Remove later
      // ------------------------------------------------------------

      Acts::Svg::Style portalStyle;
      portalStyle.strokeHighlights = {"mouseover", "mouseout"};
      portalStyle.strokeHighlightWidth = 3.;
      portalStyle.strokeHighlightColor = {255, 255, 0};

      Acts::Svg::DetectorConverter::Options detectorOptions;
      detectorOptions.volumeOptions.portalOptions.surfaceOptions.style =
          portalStyle;

      auto pDetector = Acts::Svg::DetectorConverter::convert(
          tContext, *detector, detectorOptions);
      pDetector._name = detector->name();

      for (auto& c : colors) {
        c._opacity = 0.1;
      }

      pDetector.colorize(colors);

      // As zr view
      auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
      Acts::Svg::toFile({dv_zr},
                        "DetectorVolumesInZ_" + trfStr + radStr + "_zr.svg");
    }
  }
}

BOOST_AUTO_TEST_CASE(PerfectFitPhiAttachment) {
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};
  unsigned int phiSectors = 5;
  Acts::ActsScalar phiHalfSector = M_PI / phiSectors;

  auto portalGenerator = detail::defaultPortalGenerator();
  auto navigationStateUpdator = detail::defaultPortalProvider();

  for (const auto& t : transforms) {
    std::vector<std::shared_ptr<DetectorVolume>> phiVolumes = {};
    for (unsigned int i = 0; i < phiSectors; ++i) {
      auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
          10., 100., 100., phiHalfSector, 0.);
      // Rotate in phi center
      const Acts::Vector3 colZ = t.rotation().matrix().col(2);
      Acts::Transform3 ti =
          Acts::AngleAxis3(-M_PI + 2 * i * phiHalfSector, colZ) * t;

      // create the volume
      phiVolumes.push_back(DetectorVolumeFactory::construct(
          portalGenerator, tContext, "Cylinder_phi" + std::to_string(i), ti,
          std::move(cBounds), detail::defaultPortalProvider()));
    }

    connectVolumesInPhi(tContext, phiVolumes);

    // A detector construction that should work
    DetectorVolumeFinder trialAndErrorFinder;
    trialAndErrorFinder.connect<&trialAndError>();

    auto detector =
        Detector::makeShared("DetectorInPhi", phiVolumes, trialAndErrorFinder);

    // Remove later
    // ------------------------------------------------------------
    Acts::Svg::Style portalStyle;
    portalStyle.strokeHighlights = {"mouseover", "mouseout"};
    portalStyle.strokeHighlightWidth = 3.;
    portalStyle.strokeHighlightColor = {255, 255, 0};

    Acts::Svg::DetectorConverter::Options detectorOptions;
    detectorOptions.volumeOptions.portalOptions.surfaceOptions.style =
        portalStyle;

    auto pDetector = Acts::Svg::DetectorConverter::convert(tContext, *detector,
                                                           detectorOptions);
    pDetector._name = detector->name();

    for (auto& c : colors) {
      c._opacity = 0.1;
    }

    pDetector.colorize(colors);

    // As xy view
    auto dv_xy = Acts::Svg::View::xy(pDetector, pDetector._name);
    Acts::Svg::toFile({dv_xy}, "DetectorVolumesInPhi_xy.svg");
  }
}

BOOST_AUTO_TEST_CASE(ProtoContainerRZ) {
  auto transform = Acts::Transform3::Identity();

  auto portalGenerator = detail::defaultPortalGenerator();
  auto navigationStateUpdator = detail::defaultPortalProvider();
  // A container in R
  std::vector<Acts::ActsScalar> radii = {25., 100., 200.};
  Acts::ActsScalar halfZ = 200;

  // An innermost Pipe
  auto bBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(0., radii[0u], halfZ);

  auto innerPipe = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "InnerPipe", transform, std::move(bBounds),
      detail::defaultPortalProvider());

  // Make a container representation out of it
  std::map<unsigned int, std::shared_ptr<Portal>> ipContainer;
  ipContainer[0u] = innerPipe->portalPtrs()[0u];
  ipContainer[1u] = innerPipe->portalPtrs()[1u];
  ipContainer[2u] = innerPipe->portalPtrs()[2u];

  // Create the r - sorted volumes
  std::vector<std::shared_ptr<DetectorVolume>> rVolumes = {};
  // Create the voluems
  for (auto [i, r] : Acts::enumerate(radii)) {
    if (i > 0) {
      auto cBounds =
          std::make_unique<Acts::CylinderVolumeBounds>(radii[i - 1u], r, halfZ);
      rVolumes.push_back(DetectorVolumeFactory::construct(
          portalGenerator, tContext, "Cylinder_r" + std::to_string(i),
          transform, std::move(cBounds), detail::defaultPortalProvider()));
    }
  }

  auto protoContainerInR = connectVolumesInR(tContext, rVolumes);

  std::vector<Acts::ActsScalar> zValues = {-200., -120, 10., 100., 200.};
  std::vector<std::shared_ptr<DetectorVolume>> zVolumes = {};
  for (auto [i, z] : Acts::enumerate(zValues)) {
    if (i > 0) {
      auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
          200., 300., 0.5 * (z - zValues[i - 1u]));
      // z center
      Acts::ActsScalar zCenter = 0.5 * (z + zValues[i - 1u]);
      Acts::Transform3 ti = transform;
      ti.pretranslate(transform.translation() +
                      zCenter * transform.rotation().matrix().col(2));

      // create the volume
      zVolumes.push_back(DetectorVolumeFactory::construct(
          portalGenerator, tContext, "Cylinder_z" + std::to_string(i), ti,
          std::move(cBounds), detail::defaultPortalProvider()));
    }
  }
  // Now call the connector
  auto protoContainerInZ = connectVolumesInZ(tContext, zVolumes);

  auto centralContainer = connectContainersInR(
      tContext, {ipContainer, protoContainerInR, protoContainerInZ});

  // Let's make two endcaps
  // Nec
  auto necBounds = std::make_unique<Acts::CylinderVolumeBounds>(0., 300., 50.);

  auto necTransform = Acts::Transform3::Identity();
  necTransform.pretranslate(Acts::Vector3(0., 0., -250));
  auto necVolume = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Nec", necTransform, std::move(necBounds),
      detail::defaultPortalProvider());

  std::map<unsigned int, std::shared_ptr<Portal>> necContainer;
  necContainer[0u] = necVolume->portalPtrs()[0u];
  necContainer[1u] = necVolume->portalPtrs()[1u];
  necContainer[2u] = necVolume->portalPtrs()[2u];

  // Pec container
  auto pecInnerBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(0., 175., 100.);

  auto pecOuterBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(175., 300., 100.);

  auto pecTransform = Acts::Transform3::Identity();
  pecTransform.pretranslate(Acts::Vector3(0., 0., 300));
  auto pecInner = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "PecInner", pecTransform,
      std::move(pecInnerBounds), detail::defaultPortalProvider());
  auto pecOuter = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "PecOuter", pecTransform,
      std::move(pecOuterBounds), detail::defaultPortalProvider());

  std::vector<std::shared_ptr<DetectorVolume>> pecVolumes = {pecInner,
                                                             pecOuter};

  auto pecContainer = connectVolumesInR(tContext, pecVolumes);

  // Make the full connection
  auto overallContainer = connectContainersInZ(
      tContext, {necContainer, centralContainer, pecContainer});

  //  Add them togeter
  std::vector<std::shared_ptr<DetectorVolume>> dVolumes;
  dVolumes.push_back(innerPipe);
  dVolumes.push_back(necVolume);
  dVolumes.insert(dVolumes.end(), rVolumes.begin(), rVolumes.end());
  dVolumes.insert(dVolumes.end(), zVolumes.begin(), zVolumes.end());
  dVolumes.push_back(pecInner);
  dVolumes.push_back(pecOuter);

  // A detector construction that should work
  DetectorVolumeFinder trialAndErrorFinder;
  trialAndErrorFinder.connect<&trialAndError>();

  auto detector = Detector::makeShared("DetectorFromProtoContainer", dVolumes,
                                       trialAndErrorFinder);

  // Remove later
  // ------------------------------------------------------------
  Acts::Svg::Style portalStyle;
  portalStyle.strokeHighlights = {"mouseover", "mouseout"};
  portalStyle.strokeHighlightWidth = 3.;
  portalStyle.strokeHighlightColor = {255, 255, 0};

  Acts::Svg::DetectorConverter::Options detectorOptions;
  detectorOptions.volumeOptions.portalOptions.surfaceOptions.style =
      portalStyle;

  auto pDetector = Acts::Svg::DetectorConverter::convert(tContext, *detector,
                                                         detectorOptions);
  pDetector._name = detector->name();

  for (auto& c : colors) {
    c._opacity = 0.1;
  }

  pDetector.colorize(colors);

  // As zr view
  auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
  Acts::Svg::toFile({dv_zr}, "DetectorFromProtoContainerInR_zr.svg");
}

BOOST_AUTO_TEST_SUITE_END()
