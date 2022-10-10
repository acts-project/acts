// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/detail/NavigationStateUpdators.hpp"
#include "Acts/Experimental/detail/PortalGenerators.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/ActSVG/DetectorVolumeSvgConverter.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <fstream>
#include <memory>
#include <vector>

Acts::GeometryContext tgContext;

auto nominal = Acts::Transform3::Identity();

auto portalGenerator = Acts::Experimental::detail::defaultPortalGenerator();

auto navigationStateUpdator =
    Acts::Experimental::detail::defaultPortalProvider();

BOOST_AUTO_TEST_SUITE(DetectorVolumeSvgConverter)

BOOST_AUTO_TEST_CASE(TubeCylindricalDetectorVolume) {
  // The volume definitions
  Acts::ActsScalar rInner = 10.;
  Acts::ActsScalar rOuter = 100.;
  Acts::ActsScalar zHalfL = 300.;

  Acts::Svg::Style portalStyle;
  portalStyle.fillColor = {255, 255, 255};
  portalStyle.fillOpacity = 0.;

  // A tube cylinder
  auto tubeCylinderBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(rInner, rOuter, zHalfL);

  auto tubeCylinderVolume =
      Acts::Experimental::DetectorVolumeFactory::construct(
          portalGenerator, tgContext, "TubeCylinderVolume", nominal,
          std::move(tubeCylinderBounds), navigationStateUpdator);

  Acts::Svg::DetectorVolumeConverter::Options volumeOptions;
  volumeOptions.portalOptions.volumeIndices[tubeCylinderVolume.get()] = 0u;

  auto pVolume = Acts::Svg::DetectorVolumeConverter::convert(
      tgContext, *tubeCylinderVolume, volumeOptions);
  pVolume._name = tubeCylinderVolume->name();

  // Colorize in red
  actsvg::style::color red({{255, 0, 0}});
  red._opacity = 0.1;
  std::vector<actsvg::style::color> colors = {red};
  pVolume.colorize(colors);

  // As sheet
  auto pv = Acts::Svg::View::zr(pVolume, pVolume._name);
  Acts::Svg::toFile({pv}, pVolume._name + "_zr.svg");
}

BOOST_AUTO_TEST_CASE(TubeSectorCylindricalDetectorVolume) {
  // The volume definitions
  Acts::ActsScalar rInner = 10.;
  Acts::ActsScalar rOuter = 100.;
  Acts::ActsScalar zHalfL = 300.;
  Acts::ActsScalar phiSector = 0.25 * M_PI;
  std::vector<Acts::ActsScalar> avgPhi = {0., 0.75};
  std::vector<std::string> avgPhiTag = {"zero", "nonzero"};

  Acts::Svg::Style portalStyle;
  portalStyle.fillColor = {255, 255, 255};
  portalStyle.fillOpacity = 0.;

  for (auto [iphi, aphi] : Acts::enumerate(avgPhi)) {
    // A tube cylinder
    auto sectorCylinderBounds = std::make_unique<Acts::CylinderVolumeBounds>(
        rInner, rOuter, zHalfL, phiSector, aphi);

    auto sectorCylinderVolume =
        Acts::Experimental::DetectorVolumeFactory::construct(
            portalGenerator, tgContext, "SectoralCylinderVolume", nominal,
            std::move(sectorCylinderBounds), navigationStateUpdator);

    Acts::Svg::DetectorVolumeConverter::Options volumeOptions;
    volumeOptions.portalOptions.volumeIndices[sectorCylinderVolume.get()] = 0u;

    auto pVolume = Acts::Svg::DetectorVolumeConverter::convert(
        tgContext, *sectorCylinderVolume, volumeOptions);

    // Colorize in blue
    actsvg::style::color blue({{0, 0, 255}});
    blue._opacity = 0.1;
    std::vector<actsvg::style::color> colors = {blue};
    pVolume.colorize(colors);
  }
}

BOOST_AUTO_TEST_SUITE_END()
