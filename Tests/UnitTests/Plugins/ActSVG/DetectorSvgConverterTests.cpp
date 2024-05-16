// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/CylindricalContainerBuilder.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/interface/IDetectorComponentBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Plugins/ActSVG/DetectorSvgConverter.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <fstream>
#include <memory>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;

GeometryContext tContext;

/// @brief A mockup volume builder, it generates volumes with
/// a single surface filled in in order to use the CylindricalContainerBuilder
/// infrastructure.
template <typename surface_type, typename surface_bounds_type>
class CylindricalVolumeBuilder : public IDetectorComponentBuilder {
 public:
  CylindricalVolumeBuilder(const Transform3& transform,
                           const CylinderVolumeBounds& vBounds,
                           const surface_bounds_type& sBounds,
                           const std::string& vName)
      : IDetectorComponentBuilder(),
        m_transform(transform),
        m_volumeBounds(vBounds),
        m_surfaceBounds(sBounds),
        m_name(vName) {}

  DetectorComponent construct(
      [[maybe_unused]] const GeometryContext& gctx) const final {
    // The outgoing root volumes
    std::vector<std::shared_ptr<DetectorVolume>> rootVolumes;

    // Ingredients
    auto surface = Surface::makeShared<surface_type>(
        (m_transform), std::make_shared<surface_bounds_type>(m_surfaceBounds));

    auto bounds = std::make_unique<CylinderVolumeBounds>(m_volumeBounds);
    auto portalGenerator = defaultPortalGenerator();
    auto volume = DetectorVolumeFactory::construct(
        portalGenerator, tContext, m_name, m_transform, std::move(bounds),
        {surface}, {}, tryNoVolumes(), tryAllPortalsAndSurfaces());

    // Add to the roots
    rootVolumes.push_back(volume);

    DetectorComponent::PortalContainer dContainer;
    for (auto [ip, p] : enumerate(volume->portalPtrs())) {
      dContainer[ip] = p;
    }
    return DetectorComponent{
        {volume},
        dContainer,
        RootDetectorVolumes{rootVolumes, tryRootVolumes()}};
  }

 private:
  Transform3 m_transform;
  CylinderVolumeBounds m_volumeBounds;
  surface_bounds_type m_surfaceBounds;
  std::string m_name;
};

BOOST_AUTO_TEST_SUITE(ActSvg)

BOOST_AUTO_TEST_CASE(CylindricalDetector) {
  // Declare a barrel sub builder
  auto beampipe = std::make_shared<
      CylindricalVolumeBuilder<CylinderSurface, CylinderBounds>>(
      Transform3::Identity(), CylinderVolumeBounds(0., 50., 400.),
      CylinderBounds(25., 380.), "BeamPipe");

  // Declare a negative disc builder
  Transform3 negZ = Transform3::Identity();
  negZ.pretranslate(Vector3(0., 0., -300.));
  auto endcapN =
      std::make_shared<CylindricalVolumeBuilder<DiscSurface, RadialBounds>>(
          negZ, CylinderVolumeBounds(50., 140., 100.), RadialBounds(60., 120.),
          "NegativeEndcap");

  // Declare a barrel sub builder
  auto barrel0 = std::make_shared<
      CylindricalVolumeBuilder<CylinderSurface, CylinderBounds>>(
      Transform3::Identity(), CylinderVolumeBounds(50., 80., 200.),
      CylinderBounds(65., 180.), "Barrel0");

  // Declare a barrel sub builder
  auto barrel1 = std::make_shared<
      CylindricalVolumeBuilder<CylinderSurface, CylinderBounds>>(
      Transform3::Identity(), CylinderVolumeBounds(80., 110., 200.),
      CylinderBounds(95., 180.), "Barrel1");

  // Declare a barrel sub builder
  auto barrel2 = std::make_shared<
      CylindricalVolumeBuilder<CylinderSurface, CylinderBounds>>(
      Transform3::Identity(), CylinderVolumeBounds(110., 140., 200.),
      CylinderBounds(125., 180.), "Barrel2");

  // Create the barrel container builder
  CylindricalContainerBuilder::Config barrelRCfg;
  barrelRCfg.builders = {barrel0, barrel1, barrel2};
  barrelRCfg.binning = {binR};

  auto barrel = std::make_shared<CylindricalContainerBuilder>(
      barrelRCfg, getDefaultLogger("BarrelBuilderR", Logging::VERBOSE));

  Transform3 posZ = Transform3::Identity();
  posZ.pretranslate(Vector3(0., 0., 300.));
  auto endcapP =
      std::make_shared<CylindricalVolumeBuilder<DiscSurface, RadialBounds>>(
          posZ, CylinderVolumeBounds(50., 140., 100.), RadialBounds(60., 120.),
          "PositiveEndcap");

  // Create the barrel container builder
  CylindricalContainerBuilder::Config barrelEndcapCfg;
  barrelEndcapCfg.builders = {endcapN, barrel, endcapP};
  barrelEndcapCfg.binning = {binZ};

  auto barrelEndcap = std::make_shared<CylindricalContainerBuilder>(
      barrelEndcapCfg,
      getDefaultLogger("BarrelEndcapBuilder", Logging::VERBOSE));

  // Create the barrel container builder
  CylindricalContainerBuilder::Config detectorCfg;
  detectorCfg.builders = {beampipe, barrelEndcap};
  detectorCfg.binning = {binR};

  auto containerBuilder = std::make_shared<CylindricalContainerBuilder>(
      detectorCfg, getDefaultLogger("DetectorBuilder", Logging::VERBOSE));

  // Detector builder
  auto gigConfig = GeometryIdGenerator::Config();
  auto gig = std::make_shared<GeometryIdGenerator>(gigConfig);

  Acts::Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxiliary = "*** Test : Cylindrical Detector ***";
  dCfg.name = "CylindricalDetector";
  dCfg.builder = containerBuilder;
  dCfg.geoIdGenerator = gig;

  auto detector = DetectorBuilder(dCfg).construct(tContext);

  Acts::Svg::DetectorConverter::Options detectorOptions;
  auto pDetector = Acts::Svg::DetectorConverter::convert(tContext, *detector,
                                                         detectorOptions);
  pDetector._name = detector->name();

  // Colorize in blue
  actsvg::style::color gray({{155, 155, 155}});
  actsvg::style::color pink({{255, 153, 255}});
  actsvg::style::color brown({{153, 102, 51}});
  actsvg::style::color red({{255, 0, 0}});
  actsvg::style::color green({{0, 255, 0}});
  actsvg::style::color blue({{0, 0, 255}});
  std::vector<actsvg::style::color> colors = {gray, pink,  brown,
                                              red,  green, blue};
  for (auto& c : colors) {
    c._opacity = 0.1;
  }

  pDetector.colorize(colors);

  // As sheet
  auto dv_zr = Acts::Svg::View::zr(pDetector, pDetector._name);
  Acts::Svg::toFile({dv_zr}, pDetector._name + "_zr.svg");
}

BOOST_AUTO_TEST_SUITE_END()
