// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/detail/CylindricalDetectorHelper.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "ActsPlugins/ActSVG/PortalSvgConverter.hpp"
#include "ActsPlugins/ActSVG/SvgUtils.hpp"

#include <algorithm>
#include <fstream>
#include <memory>
#include <ranges>
#include <vector>

using namespace Acts;
using namespace ActsPlugins;

GeometryContext tContext;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(ActSvgSuite)

BOOST_AUTO_TEST_CASE(CylinderPortalsSvg) {
  //  Some style parameteers
  Svg::Style portalStyle;
  portalStyle.fillColor = {255, 255, 255};
  portalStyle.fillOpacity = 0.;
  portalStyle.highlightColor = {255, 255, 255};
  portalStyle.highlights = {};
  portalStyle.strokeColor = {25, 25, 25};
  portalStyle.strokeWidth = 0.5;
  portalStyle.quarterSegments = 72u;

  double rInner = 10.;
  double rOuter = 100.;
  double zHalfL = 200.;

  Transform3 nominal = Transform3::Identity();

  auto cylinderBounds =
      std::make_unique<CylinderVolumeBounds>(rInner, rOuter, zHalfL);

  auto portalGenerator = Experimental::defaultPortalGenerator();

  auto cylinderVolume = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "CylinderVolume", nominal,
      std::move(cylinderBounds), Experimental::tryAllPortals());

  Svg::PortalConverter::Options portalOptions;
  portalOptions.volumeIndices[cylinderVolume.get()] = 0;
  Svg::SurfaceConverter::Options surfaceOptions;
  surfaceOptions.style = portalStyle;
  portalOptions.surfaceOptions = surfaceOptions;

  std::vector<Svg::ProtoPortal> protoPortals;

  std::ranges::for_each(cylinderVolume->portals(), [&](const auto& p) {
    protoPortals.push_back(
        Svg::PortalConverter::convert(tContext, *p, portalOptions));
  });

  // rz view
  std::vector<actsvg::svg::object> zrPortals;
  std::size_t ip = 0;
  std::ranges::for_each(protoPortals, [&](const auto& p) {
    zrPortals.push_back(Svg::View::zr(p, "Portal_zr" + std::to_string(ip++)));
  });
  Svg::toFile(zrPortals, "SimpleCylinderPortals_ZR.svg");

  // xy view
  std::vector<actsvg::svg::object> xyPortals;
  ip = 0;
  std::ranges::for_each(protoPortals, [&](const auto& p) {
    xyPortals.push_back(Svg::View::xy(p, "Portal_xy" + std::to_string(ip++)));
  });
  Svg::toFile(xyPortals, "SimpleCylinderPortals_XY.svg");
}

BOOST_AUTO_TEST_CASE(CylinderContainerPortalsSvg) {
  //  Some style parameteers
  Svg::Style portalStyle;
  portalStyle.fillColor = {255, 255, 255};
  portalStyle.fillOpacity = 0.;
  portalStyle.highlightColor = {255, 255, 255};
  portalStyle.highlights = {};
  portalStyle.strokeColor = {25, 25, 25};
  portalStyle.strokeWidth = 0.5;
  portalStyle.quarterSegments = 72u;

  double rInner = 10.;
  double rMiddle = 100.;
  double rOuter = 300.;
  double zHalfL = 200.;

  Transform3 nominal = Transform3::Identity();

  auto portalGenerator = Experimental::defaultPortalGenerator();

  auto cylinderBoundsI =
      std::make_unique<CylinderVolumeBounds>(rInner, rMiddle, zHalfL);

  auto cylinderBoundsO =
      std::make_unique<CylinderVolumeBounds>(rMiddle, rOuter, zHalfL);

  auto cylinderVolumeI = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "CylinderVolumeI", nominal,
      std::move(cylinderBoundsI), Experimental::tryAllPortals());

  auto cylinderVolumeO = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "CylinderVolumeO", nominal,
      std::move(cylinderBoundsO), Experimental::tryAllPortals());

  std::vector<std::shared_ptr<Experimental::DetectorVolume>> rVolumes = {
      cylinderVolumeI, cylinderVolumeO};

  Experimental::detail::CylindricalDetectorHelper::connectInR(tContext,
                                                              rVolumes, {});

  std::set<const Experimental::Portal*> portals;
  for (auto& v : rVolumes) {
    for (auto& p : v->portals()) {
      portals.insert(p);
    }
  }
  std::vector<Svg::ProtoPortal> protoPortals;

  Svg::PortalConverter::Options portalOptions;
  portalOptions.volumeIndices[cylinderVolumeI.get()] = 0;
  portalOptions.volumeIndices[cylinderVolumeO.get()] = 1;

  Svg::SurfaceConverter::Options surfaceOptions;
  surfaceOptions.style = portalStyle;
  portalOptions.surfaceOptions = surfaceOptions;

  for (const auto& portal : portals) {
    protoPortals.push_back(
        Svg::PortalConverter::convert(tContext, *portal, portalOptions));
  }
  // rz view
  std::vector<actsvg::svg::object> zrPortals;
  std::size_t ip = 0;
  std::ranges::for_each(protoPortals, [&](const auto& p) {
    zrPortals.push_back(Svg::View::zr(p, "Portal_zr" + std::to_string(ip++)));
  });
  Svg::toFile(zrPortals, "CylinderContainerPortals_ZR.svg");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
