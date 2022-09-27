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
std::vector<actsvg::style::color> colors = {red, green, blue, black, olive};

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

BOOST_AUTO_TEST_CASE(PerfectFitAttachmentR) {
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};
  std::vector<Acts::ActsScalar> radii = {0., 10., 100., 200.};
  Acts::ActsScalar halfZ = 100;

  auto portalGenerator = detail::defaultPortalGenerator();
  auto navigationStateUpdator = detail::defaultPortalProvider();

  // This should work for generic transforms
  for (auto [it, tf] : Acts::enumerate(transforms)) {
    std::vector<std::shared_ptr<DetectorVolume>> rVolumes = {};
    // Create the voluems
    for (auto [i, r] : Acts::enumerate(radii)) {
      if (i > 0) {
        auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
            radii[i - 1u], r, halfZ);
        rVolumes.push_back(DetectorVolume::makeShared(
            "Cylinder_r" + std::to_string(i), tContext, tf, std::move(cBounds),
            portalGenerator, navigationStateUpdator));
      }
    }
    // Now call the connector
    CylindricalDetectorHelperOptions cdhOptions;
    connectCylindricalVolumes(tContext, rVolumes, cdhOptions);

    // A detector construction that should work
    DetectorVolumeFinder trialAndErrorFinder;
    trialAndErrorFinder.connect<&trialAndError>();

    auto detector =
        Detector::makeShared("DetectorInR", rVolumes, trialAndErrorFinder);

    // Remove later ------------------------------------------------------------
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

    // As sheet
    auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
    Acts::Svg::toFile({dv_zr}, "CDH_rBinned _" + pDetector._name + "_zr.svg");
  }
}

BOOST_AUTO_TEST_CASE(PerfectFitAttachmentZ) {
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};
  std::vector<std::array<Acts::ActsScalar, 2>> radii = {{0., 100.}};
  std::vector<Acts::ActsScalar> zValues = {-100., -20, 10., 100., 200.};

  auto portalGenerator = detail::defaultPortalGenerator();
  auto navigationStateUpdator = detail::defaultPortalProvider();

  for (const auto& t : transforms) {
    for (const auto& r : radii) {
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
          zVolumes.push_back(DetectorVolume::makeShared(
              "Cylinder_z" + std::to_string(i), tContext, ti,
              std::move(cBounds), portalGenerator, navigationStateUpdator));
        }
      }
      // Now call the connector
      CylindricalDetectorHelperOptions cdhOptions;
      connectCylindricalVolumes(tContext, zVolumes, cdhOptions);

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

      // As sheet
      auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
      Acts::Svg::toFile({dv_zr}, "CDH_rBinned _" + pDetector._name + "_zr.svg");
    }
  }
}

BOOST_AUTO_TEST_CASE(PerfectFitPhiAttachment) {
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};
  unsigned int phiSectors = 12;
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
      phiVolumes.push_back(DetectorVolume::makeShared(
          "Cylinder_phi" + std::to_string(i), tContext, ti, std::move(cBounds),
          portalGenerator, navigationStateUpdator));
    }
    // Now call the connector
    CylindricalDetectorHelperOptions cdhOptions;
    connectCylindricalVolumes(tContext, phiVolumes, cdhOptions);
  }
}

BOOST_AUTO_TEST_SUITE_END()
