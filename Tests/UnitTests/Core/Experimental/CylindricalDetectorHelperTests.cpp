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
#include "Acts/Experimental/detail/DetectorVolumeFinders.hpp"
#include "Acts/Experimental/detail/NavigationStateUpdators.hpp"
#include "Acts/Experimental/detail/PortalGenerators.hpp"
#include "Acts/Experimental/detail/SurfaceGridGenerator.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <exception>
#include <memory>

// Remove later - debug only
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Plugins/ActSVG/DetectorSvgConverter.hpp"
#include "Acts/Plugins/Json/DetectorJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

// Colorize in blue
actsvg::style::color red =
    actsvg::style::color({{255, 0, 0}, 1., {}, {0, 0, 255}});
actsvg::style::color green =
    actsvg::style::color({{0, 255, 0}, 1., {}, {255, 0, 0}});
actsvg::style::color blue =
    actsvg::style::color({{0, 0, 255}, 1., {}, {255, 165, 0}});
actsvg::style::color black = actsvg::style::color({{0, 0, 0}});
actsvg::style::color olive = actsvg::style::color({{128, 128, 0}});
actsvg::style::color orange = actsvg::style::color({{255, 165, 0}});
actsvg::style::color cyan = actsvg::style::color({{0, 255, 255}});
actsvg::style::color magenta = actsvg::style::color({{255, 0, 255}});
actsvg::style::color brown = actsvg::style::color({{153, 102, 51}});
actsvg::style::color purple = actsvg::style::color({{153, 0, 204}});
std::vector<actsvg::style::color> colors = {
    red, green, blue, black, olive, orange, cyan, magenta, brown, purple};

using namespace Acts::Experimental;

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(ConnectPerfectlyFitInR) {
  // Test with different transforms
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};

  // Test with different opening angles
  std::vector<Acts::ActsScalar> openings = {M_PI};

  std::vector<Acts::ActsScalar> radii = {0., 10., 100., 200.};
  Acts::ActsScalar halfZ = 100.;

  auto portalGenerator = detail::defaultPortalGenerator();
  auto navigationStateUpdator = detail::allPortals();

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
              std::move(cBounds), detail::allPortals()));
        }
      }

      connectVolumesInR(tContext, rVolumes);

      // A detector construction that should work
      auto detector = Detector::makeShared("DetectorInR", rVolumes,
                                           detail::tryAllVolumes());

      // Make a rzphi grid
      const auto& volumes = detector->volumes();
      auto boundaries = rzphiBoundaries(tContext, volumes);
      const auto& rBoundaries = boundaries[0u];
      const auto& zBoundaries = boundaries[1u];

      // Check the radii
      std::vector<Acts::ActsScalar> zvalues = {-halfZ, halfZ};
      BOOST_CHECK(radii == rBoundaries);
      BOOST_CHECK(zvalues == zBoundaries);

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
  auto navigationStateUpdator = detail::allPortals();

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
              std::move(cBounds), detail::allPortals()));
        }
      }
      // Now call the connector
      connectVolumesInZ(tContext, zVolumes);

      auto detector = Detector::makeShared("DetectorInZ", zVolumes,
                                           detail::tryAllVolumes());

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
  auto navigationStateUpdator = detail::allPortals();

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
          std::move(cBounds), detail::allPortals()));
    }

    connectVolumesInPhi(tContext, phiVolumes);
    auto detector = Detector::makeShared("DetectorInPhi", phiVolumes,
                                         detail::tryAllVolumes());

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

  std::vector<Acts::ActsScalar> innerMostRadii = {0., 2.};

  for (auto [ir, imr] : Acts::enumerate(innerMostRadii)) {
    auto portalGenerator = detail::defaultPortalGenerator();
    auto navigationStateUpdator = detail::allPortals();
    // A container in R
    std::vector<Acts::ActsScalar> radii = {25., 100., 200.};
    Acts::ActsScalar halfZ = 200;

    // An innermost Pipe
    auto bBounds =
        std::make_unique<Acts::CylinderVolumeBounds>(imr, radii[0u], halfZ);

    auto innerPipe = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "InnerPipe", transform, std::move(bBounds),
        detail::allPortals());

    // Make a container representation out of it
    std::map<unsigned int, std::shared_ptr<Portal>> ipContainer;
    for (auto [ip, p] : Acts::enumerate(innerPipe->portalPtrs())) {
      ipContainer[ip] = p;
    }

    // Create the r - sorted volumes
    std::vector<std::shared_ptr<DetectorVolume>> rVolumes = {};
    // Create the voluems
    for (auto [i, r] : Acts::enumerate(radii)) {
      if (i > 0) {
        auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
            radii[i - 1u], r, halfZ);
        rVolumes.push_back(DetectorVolumeFactory::construct(
            portalGenerator, tContext, "Cylinder_r" + std::to_string(i),
            transform, std::move(cBounds), detail::allPortals()));
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
            std::move(cBounds), detail::allPortals()));
      }
    }
    // Now call the connector
    auto protoContainerInZ = connectVolumesInZ(tContext, zVolumes);

    auto centralContainer = connectContainersInR(
        tContext, {ipContainer, protoContainerInR, protoContainerInZ});

    // Let's make two endcaps
    // Nec
    auto necBounds =
        std::make_unique<Acts::CylinderVolumeBounds>(imr, 300., 50.);

    auto necTransform = Acts::Transform3::Identity();
    necTransform.pretranslate(Acts::Vector3(0., 0., -250));
    auto necVolume = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "Nec", necTransform, std::move(necBounds),
        detail::allPortals());

    std::map<unsigned int, std::shared_ptr<Portal>> necContainer;
    for (auto [ip, p] : Acts::enumerate(necVolume->portalPtrs())) {
      necContainer[ip] = p;
    }

    // Pec container
    auto pecInnerBounds =
        std::make_unique<Acts::CylinderVolumeBounds>(imr, 175., 100.);

    auto pecOuterBounds =
        std::make_unique<Acts::CylinderVolumeBounds>(175., 300., 100.);

    auto pecTransform = Acts::Transform3::Identity();
    pecTransform.pretranslate(Acts::Vector3(0., 0., 300));
    auto pecInner = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "PecInner", pecTransform,
        std::move(pecInnerBounds), detail::allPortals());
    auto pecOuter = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "PecOuter", pecTransform,
        std::move(pecOuterBounds), detail::allPortals());

    std::vector<std::shared_ptr<DetectorVolume>> pecVolumes = {pecInner,
                                                               pecOuter};
    auto pecContainer = connectVolumesInR(tContext, pecVolumes);

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

    auto detector = Detector::makeShared("DetectorFromProtoContainer", dVolumes,
                                         detail::tryAllVolumes());

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
    Acts::Svg::toFile({dv_zr}, "DetectorFromProtoContainers_" +
                                   std::to_string(ir) + "_zr.svg");
  }
}

BOOST_AUTO_TEST_CASE(ODDFrame) {
  // Standard portal generator, context, transform
  auto portalGenerator = detail::defaultPortalGenerator();
  Acts::GeometryContext tContext;
  auto transform = Acts::Transform3::Identity();

  // Volumes & Coloring
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> dVolumes =
      {};
  std::vector<actsvg::style::color> oddColors = {};

  // Beam Pipe radius
  Acts::ActsScalar bpR = 23.;
  Acts::ActsScalar bpOuterR = 28.;
  Acts::ActsScalar detectorHalfZ = 3100.;
  auto bp = Acts::Surface::makeShared<Acts::CylinderSurface>(
      Acts::Transform3::Identity(), bpR, detectorHalfZ - 10.);
  std::vector<std::shared_ptr<Acts::Surface>> bpSurfaces = {bp};
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> empty = {};

  auto bpBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(0, bpOuterR, detectorHalfZ);

  auto bpVolume = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "BeamPipe", transform, std::move(bpBounds),
      bpSurfaces, empty, detail::allPortalsAndSurfaces());

  // Make a container representation out of it
  std::map<unsigned int, std::shared_ptr<Portal>> bpContainer;
  for (auto [ip, p] : Acts::enumerate(bpVolume->portalPtrs())) {
    bpContainer[ip] = p;
  }
  // Get those volumes
  dVolumes.push_back(bpVolume);
  oddColors.push_back(orange);
  size_t cVolumes = dVolumes.size();

  // EC helper function -> will go to Cylindrical Detector helper
  auto ecContainer =
      [&](const std::vector<Acts::ActsScalar>& zValues, Acts::ActsScalar innerR,
          Acts::ActsScalar outerR,
          const std::string& nameTag) -> Acts::Experimental::ProtoContainer {
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> zVolumes;
    for (auto [i, z] : Acts::enumerate(zValues)) {
      if (i > 0) {
        auto bounds = std::make_unique<Acts::CylinderVolumeBounds>(
            innerR, outerR, 0.5 * (z - zValues[i - 1u]));
        // z center
        Acts::ActsScalar zCenter = 0.5 * (z + zValues[i - 1u]);
        Acts::Transform3 ti = Acts::Transform3::Identity();
        ti.pretranslate(transform.translation() +
                        zCenter * transform.rotation().matrix().col(2));
        ti.prerotate(transform.rotation());
        // create the volume
        zVolumes.push_back(DetectorVolumeFactory::construct(
            portalGenerator, tContext, nameTag + "_vol_" + std::to_string(i),
            ti, std::move(bounds), detail::allPortals()));
      }
    }
    // Insert the volumes
    dVolumes.insert(dVolumes.end(), zVolumes.begin(), zVolumes.end());
    return connectVolumesInZ(tContext, zVolumes);
  };

  // Barrel helper function -> will go to Cylindrical Detector helper
  auto bContainer =
      [&](const std::vector<Acts::ActsScalar>& rValues, Acts::ActsScalar minZ,
          Acts::ActsScalar maxZ,
          const std::string& nameTag) -> Acts::Experimental::ProtoContainer {
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> rVolumes;
    Acts::ActsScalar zCenter = 0.5 * (minZ + maxZ);
    Acts::ActsScalar halfZ = 0.5 * (maxZ - minZ);
    Acts::Transform3 ti = Acts::Transform3::Identity();
    ti.pretranslate(transform.translation() +
                    zCenter * transform.rotation().matrix().col(2));
    ti.prerotate(transform.rotation());

    // Create the volumes
    for (auto [i, r] : Acts::enumerate(rValues)) {
      if (i > 0) {
        auto bounds = std::make_unique<Acts::CylinderVolumeBounds>(
            rValues[i - 1u], r, halfZ);
        rVolumes.push_back(DetectorVolumeFactory::construct(
            portalGenerator, tContext, nameTag + "_vol_" + std::to_string(i),
            ti, std::move(bounds), detail::allPortals()));
      }
    }
    // Insert the volumes
    dVolumes.insert(dVolumes.end(), rVolumes.begin(), rVolumes.end());
    return connectVolumesInR(tContext, rVolumes);
  };

  // Pixel pack
  Acts::ActsScalar pixInnerR = bpOuterR;
  Acts::ActsScalar pixOuterR = 195.;
  Acts::ActsScalar pixBoundaryECB = 600.;
  std::vector<Acts::ActsScalar> ecBoundaries = {
      -detectorHalfZ, -1970, -1930,  -1540., -1500.,         -1340.,
      -1300.,         -1140, -1100., -1000., -960.,          -860.,
      -820.,          -740., -700.,  -640.,  -pixBoundaryECB};
  std::vector<Acts::ActsScalar> bBoundaries = {
      pixInnerR, 39., 65., 77., 100., 126., 160., 182., pixOuterR};

  // Get those negative endcap volumes
  auto pixNecContainer =
      ecContainer(ecBoundaries, pixInnerR, pixOuterR, "PixelNegativeEndcap");
  // Get those central volumes
  auto pixBContainer =
      bContainer(bBoundaries, -pixBoundaryECB, pixBoundaryECB, "PixelBarrel");
  // Get those negative endcap volumes
  std::for_each(ecBoundaries.begin(), ecBoundaries.end(),
                [&](auto& s) { s *= -1; });

  std::sort(ecBoundaries.begin(), ecBoundaries.end());
  auto pixPecContainer =
      ecContainer(ecBoundaries, pixInnerR, pixOuterR, "PixelPositiveEndcap");

  auto pixContainer = connectContainersInZ(
      tContext, {pixNecContainer, pixBContainer, pixPecContainer});
  for (unsigned int dv = cVolumes; dv < dVolumes.size(); ++dv) {
    oddColors.push_back(blue);
  }
  cVolumes = dVolumes.size();

  // Pixel support tube
  Acts::ActsScalar pstInnerR = pixOuterR;
  Acts::ActsScalar pstOuterR = 205.;
  auto pstBounds = std::make_unique<Acts::CylinderVolumeBounds>(
      pstInnerR, pstOuterR, detectorHalfZ);
  auto pstVolume = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "PST", transform, std::move(pstBounds),
      detail::allPortals());
  // Make a container representation out of it
  std::map<unsigned int, std::shared_ptr<Portal>> pstContainer;
  for (auto [ip, p] : Acts::enumerate(pstVolume->portalPtrs())) {
    pstContainer[ip] = p;
  }
  // Add pst
  dVolumes.push_back(pstVolume);
  oddColors.push_back(purple);
  cVolumes = dVolumes.size();

  // Short strips
  Acts::ActsScalar ssInnerR = pstOuterR;
  Acts::ActsScalar ssOuterR = 730.;
  Acts::ActsScalar sBoundaryECB = 1265.;

  ecBoundaries = {-detectorHalfZ, -2925, -2575,  -2525., -2225., -2175.,
                  -1875.,         -1825, -1575., -1525., -1325,  -sBoundaryECB};

  bBoundaries = {ssInnerR, 240., 280., 340., 380.,
                 480.,     520., 640., 680., ssOuterR};

  // Get those negative endcap volumes
  auto ssNecContainer = ecContainer(ecBoundaries, ssInnerR, ssOuterR,
                                    "ShortStripsNegativeEndcap");
  // Get those central volumes
  auto ssBContainer =
      bContainer(bBoundaries, -sBoundaryECB, sBoundaryECB, "ShortStripsBarrel");
  // Get those negative endcap volumes
  std::for_each(ecBoundaries.begin(), ecBoundaries.end(),
                [&](auto& s) { s *= -1; });

  std::sort(ecBoundaries.begin(), ecBoundaries.end());
  auto ssPecContainer = ecContainer(ecBoundaries, ssInnerR, ssOuterR,
                                    "ShortStripsPositiveEndcap");

  auto ssContainer = connectContainersInZ(
      tContext, {ssNecContainer, ssBContainer, ssPecContainer});
  for (unsigned int dv = cVolumes; dv < dVolumes.size(); ++dv) {
    oddColors.push_back(red);
  }
  cVolumes = dVolumes.size();

  // Long strips
  Acts::ActsScalar lsInnerR = ssOuterR;
  Acts::ActsScalar lsOuterR = 1080.;

  ecBoundaries = {-detectorHalfZ, -2965, -2645,  -2555., -2295., -2205.,
                  -1945.,         -1855, -1645., -1555., -1345,  -sBoundaryECB};

  bBoundaries = {lsInnerR, 790., 840., 990., 1040., lsOuterR};

  // Get those negative endcap volumes
  auto lsNecContainer =
      ecContainer(ecBoundaries, lsInnerR, lsOuterR, "LongStripsNegativeEndcap");
  // Get those central volumes
  auto lsBContainer =
      bContainer(bBoundaries, -sBoundaryECB, sBoundaryECB, "LongStripsBarrel");
  // Get those negative endcap volumes
  std::for_each(ecBoundaries.begin(), ecBoundaries.end(),
                [&](auto& s) { s *= -1; });

  std::sort(ecBoundaries.begin(), ecBoundaries.end());
  auto lsPecContainer =
      ecContainer(ecBoundaries, lsInnerR, lsOuterR, "LongStripsPositiveEndcap");

  auto lsContainer = connectContainersInZ(
      tContext, {lsNecContainer, lsBContainer, lsPecContainer});
  for (unsigned int dv = cVolumes; dv < dVolumes.size(); ++dv) {
    oddColors.push_back(green);
  }
  cVolumes = dVolumes.size();

  // Solenoid volume
  Acts::ActsScalar solenoidInnerR = lsOuterR;
  Acts::ActsScalar solenoidOuterR = 1200.;
  auto solenoidBounds = std::make_unique<Acts::CylinderVolumeBounds>(
      solenoidInnerR, solenoidOuterR, detectorHalfZ);
  auto solenoidVolume = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Solenoid", transform,
      std::move(solenoidBounds), detail::allPortals());
  // Make a container representation out of it
  std::map<unsigned int, std::shared_ptr<Portal>> solenoidContainer;
  for (auto [ip, p] : Acts::enumerate(solenoidVolume->portalPtrs())) {
    solenoidContainer[ip] = p;
  }
  // Add pst
  dVolumes.push_back(solenoidVolume);
  oddColors.push_back(brown);
  cVolumes = dVolumes.size();

  // Build the container
  auto detectorContainer = connectContainersInR(
      tContext, {bpContainer, pixContainer, pstContainer, ssContainer,
                 lsContainer, solenoidContainer});

  // Let's assign the senstive surfaces
  auto oddIn = std::ifstream("odd_surfaces_in.json");
  nlohmann::json jidd;
  oddIn >> jidd;
  auto jSurfaces = jidd["entries"];

  std::vector<std::shared_ptr<Acts::Surface>> dSurfaces;
  for (const auto& js : jSurfaces) {
    int sIndex = js["sensitive"];
    if (sIndex > 0) {
      dSurfaces.push_back(Acts::surfaceFromJson(js["value"]));
    }
  }

  // Add the pixel endplates
  Acts::Transform3 ppepT = Acts::Transform3::Identity();
  ppepT.pretranslate(Acts::Vector3(0., 0., 1950.));
  auto ppep = Acts::Surface::makeShared<Acts::DiscSurface>(
      ppepT, pixInnerR + 10., pixOuterR - 10.);
  dSurfaces.push_back(ppep);

  Acts::Transform3 npepT = Acts::Transform3::Identity();
  npepT.pretranslate(Acts::Vector3(0., 0., -1950.));
  auto npep = Acts::Surface::makeShared<Acts::DiscSurface>(
      npepT, pixInnerR + 10., pixOuterR - 10.);
  dSurfaces.push_back(npep);

  // Add the pst
  auto pst = Acts::Surface::makeShared<Acts::CylinderSurface>(
      Acts::Transform3::Identity(), 200., detectorHalfZ - 10.);
  dSurfaces.push_back(pst);

  // Add the solenoid
  auto sol = Acts::Surface::makeShared<Acts::CylinderSurface>(
      Acts::Transform3::Identity(), 1150., 2900.);
  dSurfaces.push_back(sol);

  unsigned int assigned = 0;
  std::map<std::shared_ptr<Acts::Experimental::DetectorVolume>,
           std::vector<std::shared_ptr<Acts::Surface>>>
      assignedSurfaces;
  for (auto s : dSurfaces) {
    const auto sCenter = s->binningPosition(tContext, Acts::binR);
    for (auto v : dVolumes) {
      if (v->inside(tContext, sCenter)) {
        if (assignedSurfaces.find(v) == assignedSurfaces.end()) {
          assignedSurfaces[v] = {};
        }
        assignedSurfaces[v].push_back(s);
        ++assigned;
        break;
      }
    }
  }

  for (auto v : dVolumes) {
    if (assignedSurfaces.find(v) != assignedSurfaces.end()) {
      // Get the surfaces and build a proto-layer
      auto surfaces = assignedSurfaces.find(v)->second;
      ManagedNavigationStateUpdator mUpdator;
      std::array<Acts::BinningValue, 2u> bValues = {Acts::binR, Acts::binPhi};
      std::array<Acts::ActsScalar, 2u> tolerances = {15., 0.05};
      std::array<Acts::ActsScalar, 2u> diffTolerances = {7.5, 0.045};
      bool phiCenter = false;
      if (v->name().find("Endcap") == std::string::npos) {
        bValues = {Acts::binZ, Acts::binPhi};
        phiCenter = true;
      }
      Acts::Experimental::detail::generateNavigationStateUpdator(
          tContext, mUpdator, surfaces, bValues, tolerances, diffTolerances,
          phiCenter);
      // Update the navigaiton state updator
      v->updateNavigationStateUpator(std::move(mUpdator), surfaces);
    }
  }

  auto detector =
      Detector::makeShared("ODD", dVolumes, detail::tryAllVolumes());

  // Convert to json --- standard format
  nlohmann::json jodd;
  jodd["detector"] = toJson(*detector.get());

  std::ofstream oddOut;
  oddOut.open("ODD.json");
  oddOut << jodd.dump(4);
  oddOut.close();

  // Convert to json --- detray format
  nlohmann::json jodd_detray;
  jodd_detray["detector"] = toJson(*detector.get(), Acts::GeometryContext(),
                                   true, Acts::Logging::INFO);

  std::ofstream oddOut_detray;
  oddOut_detray.open("ODD_detray.json");
  oddOut_detray << jodd_detray.dump(4);
  oddOut_detray.close();

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

  for (auto& c : oddColors) {
    c._opacity = 0.1;
  }

  pDetector.colorize(oddColors);

  // As zr view
  auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
  Acts::Svg::toFile({dv_zr}, "ODD_zr.svg");

  // Loop over valume and make layer views:
  for (const auto& volume : pDetector._volumes) {
    if (not volume._surfaces.empty()) {
      auto [module_sheet, grid_sheet] = Acts::Svg::View::layer(volume);
      Acts::Svg::toFile({module_sheet}, volume._name + "_module_sheet.svg");
      Acts::Svg::toFile({grid_sheet}, volume._name + "_grid_sheet.svg");
    }
  }
}

BOOST_AUTO_TEST_CASE(ReadODD) {
  // Let's assign the senstive surfaces
  auto oddIn = std::ifstream("ODD.json");
  nlohmann::json jodd;
  oddIn >> jodd;

  auto detector = detectorFromJson(jodd);

  std::vector<actsvg::style::color> oddColors =
      std::vector<actsvg::style::color>(detector->volumes().size(), blue);
  // No time for colorization
  for (auto& c : oddColors) {
    c._opacity = 0.1;
  }

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

  for (auto& c : oddColors) {
    c._opacity = 0.1;
  }
  pDetector.colorize(oddColors);

  // As zr view
  auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
  Acts::Svg::toFile({dv_zr}, "ODD_read_in_zr.svg");

  // Loop over valume and make layer views:
  for (const auto& volume : pDetector._volumes) {
    if (not volume._surfaces.empty()) {
      auto [module_sheet, grid_sheet] = Acts::Svg::View::layer(volume);
      Acts::Svg::toFile({grid_sheet}, volume._name + "_grid_sheet_read_in.svg");
    }
  }

  // Convert to json --- standard format
  nlohmann::json jodd_out;
  jodd_out["detector"] = toJson(*detector.get());

  std::ofstream oddOut;
  oddOut.open("ODD_in_out.json");
  oddOut << jodd_out.dump(4);
  oddOut.close();
}

BOOST_AUTO_TEST_CASE(ReadODD_detray) {
  // Let's assign the senstive surfaces
  auto oddIn = std::ifstream("ODD_detray.json");
  nlohmann::json jodd;
  oddIn >> jodd;

  auto detector = detectorFromJson(jodd, Acts::Logging::INFO);

  std::vector<actsvg::style::color> oddColors =
      std::vector<actsvg::style::color>(detector->volumes().size(), blue);
  // No time for colorization
  for (auto& c : oddColors) {
    c._opacity = 0.1;
  }

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

  for (auto& c : oddColors) {
    c._opacity = 0.1;
  }
  pDetector.colorize(oddColors);

  // As zr view
  auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
  Acts::Svg::toFile({dv_zr}, "DETRAY_read_in_ODD_zr.svg");

  // Loop over valume and make layer views:
  for (const auto& volume : pDetector._volumes) {
    if (not volume._surfaces.empty()) {
      auto [module_sheet, grid_sheet] = Acts::Svg::View::layer(volume);
      Acts::Svg::toFile({grid_sheet},
                        "DETRAY_read_in_" + volume._name + "_grid_sheet.svg");
    }
    auto v_zr = Acts::Svg::View::zr(volume, volume._name);
    Acts::Svg::toFile({v_zr}, "DETRAY_read_in_" + volume._name + ".svg");
  }
}

BOOST_AUTO_TEST_SUITE_END()
