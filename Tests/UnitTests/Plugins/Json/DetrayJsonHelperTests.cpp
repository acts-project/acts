// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/CylindricalDetectorHelper.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/detail/NavigationStateUpdators.hpp"
#include "Acts/Experimental/detail/PortalGenerators.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/DetrayJsonHelper.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <fstream>
#include <iostream>
#include <string>

// Remove later
#include "Acts/Experimental/detail/DetectorVolumeFinders.hpp"
#include "Acts/Plugins/ActSVG/DetectorSvgConverter.hpp"

using namespace Acts;

GeometryContext tContext;

auto portalGenerator = Experimental::detail::defaultPortalGenerator();

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

namespace {

/// Generate volumes in z
///
/// @param rMinMax the radii
/// @param boundaries the boundaris in Z
std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumesInZ(
    const std::string& baseName, std::array<ActsScalar, 2u> rMinMax,
    const std::vector<ActsScalar>& boundaries) {
  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes;
  for (auto [i, z] : Acts::enumerate(boundaries)) {
    if (i > 0) {
      auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
          rMinMax[0u], rMinMax[1u], 0.5 * (z - boundaries[i - 1u]));
      // z center
      Acts::ActsScalar zCenter = 0.5 * (z + boundaries[i - 1u]);
      Acts::Transform3 ti = Acts::Transform3::Identity();
      ti.pretranslate(zCenter * ti.rotation().matrix().col(2));
      // create the volume
      volumes.push_back(Experimental::DetectorVolumeFactory::construct(
          portalGenerator, tContext, baseName + std::to_string(i), ti,
          std::move(cBounds), Experimental::detail::allPortals()));
    }
  }
  return volumes;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(DetrayJsonHelper)

BOOST_AUTO_TEST_CASE(BasicClippingTest) {
  std::vector<ActsScalar> boundaries = {-100., 200., 300., 500.};
  auto volumes = volumesInZ("Volumes_", {0., 100.}, boundaries);

  std::vector<const Experimental::DetectorVolume*> cVolumes =
      unpack_shared_const_vector(volumes);

  // test - range bigger than volume range
  std::array<ActsScalar, 2u> clipRange = {-200, 700.};
  auto [b0, v0] = Experimental::DetrayJsonHelper::clip(
      boundaries, cVolumes, clipRange, Acts::Logging::VERBOSE);
  BOOST_CHECK(b0 == boundaries);
  BOOST_CHECK(v0 == cVolumes);
  // test - range equal to volume range
  clipRange = {-100, 500.};
  auto [b1, v1] = Experimental::DetrayJsonHelper::clip(
      boundaries, cVolumes, clipRange, Acts::Logging::VERBOSE);
  BOOST_CHECK(b1 == boundaries);
  BOOST_CHECK(v1 == cVolumes);
  // test - range smaller, but no volume dropped at low
  clipRange = {-50, 500.};
  auto [b2, v2] = Experimental::DetrayJsonHelper::clip(
      boundaries, cVolumes, clipRange, Acts::Logging::VERBOSE);
  BOOST_CHECK(b2.size() == boundaries.size());
  BOOST_CHECK(b2[0u] == clipRange[0u]);
  BOOST_CHECK(v2 == cVolumes);
  // test - removing the first volume
  clipRange = {200., 500.};
  auto [b3, v3] = Experimental::DetrayJsonHelper::clip(
      boundaries, cVolumes, clipRange, Acts::Logging::VERBOSE);
  BOOST_CHECK(b3.size() == (boundaries.size() - 1u));
  BOOST_CHECK(b3[0u] == clipRange[0u]);
  BOOST_CHECK(v3 != cVolumes);
  BOOST_CHECK(std::find(v3.begin(), v3.end(), cVolumes[0u]) == v3.end());
  // test - removing the frist two volumes with setting new range of third
  clipRange = {350., 500.};
  auto [b4, v4] = Experimental::DetrayJsonHelper::clip(
      boundaries, cVolumes, clipRange, Acts::Logging::VERBOSE);
  BOOST_CHECK(b4.size() == (boundaries.size() - 2u));
  BOOST_CHECK(b4[0u] == clipRange[0u]);
  BOOST_CHECK(v4 != cVolumes);
  BOOST_CHECK(std::find(v4.begin(), v4.end(), cVolumes[0u]) == v4.end());
  BOOST_CHECK(std::find(v4.begin(), v4.end(), cVolumes[1u]) == v4.end());
  // test - removing last volume
  clipRange = {-100., 300.};
  auto [b5, v5] = Experimental::DetrayJsonHelper::clip(
      boundaries, cVolumes, clipRange, Acts::Logging::VERBOSE);
  BOOST_CHECK(b5.size() == (boundaries.size() - 1u));
  BOOST_CHECK(b5.back() == clipRange[1u]);
  BOOST_CHECK(v5 != cVolumes);
  BOOST_CHECK(std::find(v5.begin(), v5.end(), cVolumes.back()) == v5.end());
  // test - removing last two volumes and setting new end range
  clipRange = {-100., 100.};
  auto [b6, v6] = Experimental::DetrayJsonHelper::clip(
      boundaries, cVolumes, clipRange, Acts::Logging::VERBOSE);
  BOOST_CHECK(b6.size() == (boundaries.size() - 2u));
  BOOST_CHECK(b6.back() == clipRange[1u]);
  BOOST_CHECK(v6 != cVolumes);
  BOOST_CHECK(std::find(v6.begin(), v6.end(), cVolumes[cVolumes.size() - 2u]) ==
              v6.end());
  BOOST_CHECK(std::find(v6.begin(), v6.end(), cVolumes.back()) == v6.end());
  // test - remove from both sides
  clipRange = {210., 290.};
  auto [b7, v7] = Experimental::DetrayJsonHelper::clip(
      boundaries, cVolumes, clipRange, Acts::Logging::VERBOSE);
  BOOST_CHECK(b7.size() == (boundaries.size() - 2u));
  BOOST_CHECK(b7[0u] == clipRange[0u]);
  BOOST_CHECK(b7.back() == clipRange[1u]);
  BOOST_CHECK(v7 != cVolumes);
  BOOST_CHECK(std::find(v7.begin(), v7.end(), cVolumes[0u]) == v7.end());
  BOOST_CHECK(std::find(v7.begin(), v7.end(), cVolumes.back()) == v7.end());
}

BOOST_AUTO_TEST_CASE(PortalClippingAndPickingTestZ) {

  std::cout << " ------------- start from here ------------------ " << std::endl;  
  // Create volumes in z
  std::vector<ActsScalar> innerBoundaries = {-100., 200., 300., 500.};
  std::vector<ActsScalar> outerBoundaries = {-100., 400., 500.};
  auto innerVolumes = volumesInZ("InnerVolumes_", {10., 100.}, innerBoundaries);
  auto outerVolumes =
      volumesInZ("OuterVolumes_", {100., 200.}, outerBoundaries);

  // First create container
  auto innerContainer = Experimental::connectVolumesInZ(tContext, innerVolumes);
  auto outerContainer = Experimental::connectVolumesInZ(tContext, outerVolumes);
  auto fullContainer = Experimental::connectContainersInR(
      tContext, {innerContainer, outerContainer});

  // Draw it -------------------------

  Acts::Svg::DetectorConverter::Options detectorOptions;

  auto dVolumes = innerVolumes;
  dVolumes.insert(dVolumes.end(), outerVolumes.begin(), outerVolumes.end());

  // Make a detector
  auto detector = Experimental::Detector::makeShared("SomeVolumes", dVolumes,
                                                      Experimental::detail::tryAllVolumes());

  auto pDetector =
      Svg::DetectorConverter::convert(tContext, *detector, detectorOptions);
  pDetector._name = detector->name();

  for (auto& c : colors) {
    c._opacity = 0.1;
  }
  pDetector.colorize(colors);

  // As zr view
  auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
  Svg::toFile({dv_zr}, "PRECLIPPING_inZ.svg");

  // ---------------------------------------

  Svg::PortalConverter::Options pcOptions;
  for (auto [iv, dv] : enumerate(dVolumes)) {
    pcOptions.volumeIndices[dv.get()] = iv;
  }

  // Now clip them all
  std::vector<actsvg::svg::object> portalsZR;
  for (auto& dv : dVolumes) {
    auto clippedPortals = Experimental::DetrayJsonHelper::clipAndPickPortals(
        tContext, *dv.get(), Logging::VERBOSE);
    for (auto [ip, cp] : enumerate(clippedPortals)) {
      auto pPortal = Svg::PortalConverter::convert(tContext, *cp, pcOptions);
      portalsZR.push_back(
          Svg::View::zr(pPortal, dv->name() + "_portal_" + std::to_string(ip)));
    }
  }
  // As zr view
  Svg::toFile(portalsZR, "POSTCLIPPING_inZ.svg");
}

BOOST_AUTO_TEST_SUITE_END()