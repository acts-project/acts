// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
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
#include "Acts/Plugins/Detray/DetrayGeometryConverter.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <vector>

#include <detray/io/frontend/payloads.hpp>

using namespace Acts;
using namespace Acts::Experimental;

GeometryContext tContext;

auto logger =
    Acts::getDefaultLogger("DetrayGeometryConverterTests", Acts::Logging::INFO);

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

BOOST_AUTO_TEST_SUITE(DetrayConversion)

BOOST_AUTO_TEST_CASE(DetrayTransformConversion) {
  auto transform = Transform3::Identity();
  transform.pretranslate(Vector3(1., 2., 3.));
  transform.rotate(Eigen::AngleAxisd(M_PI / 2., Vector3::UnitZ()));

  detray::io::transform_payload payload =
      DetrayGeometryConverter::convertTransform(transform);
  // Transform is correctly translated
  CHECK_CLOSE_ABS(payload.tr[0u], 1.,
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(payload.tr[1u], 2.,
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(payload.tr[2u], 3.,
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  // Rotation is correctly translated
  const auto rotation = transform.rotation();
  CHECK_CLOSE_ABS(payload.rot[0u], rotation(0, 0),
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[1u], rotation(0, 1),
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[2u], rotation(0, 2),
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[3u], rotation(1, 0),
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[4u], rotation(1, 1),
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[5u], rotation(1, 2),
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[6u], rotation(2, 0),
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[7u], rotation(2, 1),
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(payload.rot[8u], rotation(2, 2),
                  std::numeric_limits<Acts::ActsScalar>::epsilon());
}

BOOST_AUTO_TEST_CASE(DetrayMaskConversion) {
  // Placeholder, masks are tested for the moment through the DetrayJsonHelper,
  // this code will move over here when the "toJsonDetray" will become
  // deprecated
}

BOOST_AUTO_TEST_CASE(DetraySurfaceConversion) {
  // Translate a cylinder
  auto cylinderSurface = Acts::Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(20., 100.));

  auto sgID = Acts::GeometryIdentifier().setSensitive(1);
  cylinderSurface->assignGeometryId(sgID);

  detray::io::surface_payload payload = DetrayGeometryConverter::convertSurface(
      tContext, *cylinderSurface, false);

  // Check the payload
  BOOST_CHECK(!payload.index_in_coll.has_value());
  BOOST_CHECK(payload.mask.shape == detray::io::shape_id::cylinder2);
  BOOST_CHECK_EQUAL(payload.source, sgID.value());
  BOOST_CHECK(payload.type == detray::surface_id::e_sensitive);
}

BOOST_AUTO_TEST_CASE(DetrayVolumeConversion) {
  auto beampipe = std::make_shared<
      CylindricalVolumeBuilder<CylinderSurface, CylinderBounds>>(
      Transform3::Identity(), CylinderVolumeBounds(0., 50., 400.),
      CylinderBounds(25., 380.), "BeamPipe");

  auto [volumes, portals, rootVolumes] = beampipe->construct(tContext);
  auto volume = volumes.front();

  DetrayConversionUtils::GeometryIdCache geoIdCache;

  detray::io::volume_payload payload = DetrayGeometryConverter::convertVolume(
      geoIdCache, tContext, *volume, {volume.get()}, *logger);

  // Check the volume payload
  BOOST_CHECK(payload.name == "BeamPipe");
  BOOST_CHECK(payload.type == detray::volume_id::e_cylinder);
  // 3 portals and 1 surface contained
  BOOST_CHECK_EQUAL(payload.surfaces.size(), 4u);
  // Let's count
  std::size_t nPortals = 0;
  for (auto& s : payload.surfaces) {
    if (s.type == detray::surface_id::e_portal) {
      nPortals++;
    }
  }
  BOOST_CHECK_EQUAL(nPortals, 3u);
  // No acceleration structure for the moment
  BOOST_CHECK(!payload.acc_links.has_value());
}

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
  barrelRCfg.binning = {BinningValue::binR};

  auto barrel = std::make_shared<CylindricalContainerBuilder>(
      barrelRCfg, getDefaultLogger("BarrelBuilderR", Logging::INFO));

  Transform3 posZ = Transform3::Identity();
  posZ.pretranslate(Vector3(0., 0., 300.));
  auto endcapP =
      std::make_shared<CylindricalVolumeBuilder<DiscSurface, RadialBounds>>(
          posZ, CylinderVolumeBounds(50., 140., 100.), RadialBounds(60., 120.),
          "PositiveEndcap");

  // Create the barrel container builder
  CylindricalContainerBuilder::Config barrelEndcapCfg;
  barrelEndcapCfg.builders = {endcapN, barrel, endcapP};
  barrelEndcapCfg.binning = {BinningValue::binZ};

  auto barrelEndcap = std::make_shared<CylindricalContainerBuilder>(
      barrelEndcapCfg, getDefaultLogger("BarrelEndcapBuilder", Logging::INFO));

  // Create the barrel container builder
  CylindricalContainerBuilder::Config detectorCfg;
  detectorCfg.builders = {beampipe, barrelEndcap};
  detectorCfg.binning = {BinningValue::binR};

  auto containerBuilder = std::make_shared<CylindricalContainerBuilder>(
      detectorCfg, getDefaultLogger("DetectorBuilder", Logging::INFO));

  // Detector builder
  auto gigConfig = GeometryIdGenerator::Config();
  auto gig = std::make_shared<GeometryIdGenerator>(gigConfig);

  Acts::Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxiliary = "*** Test : Cylindrical Detector ***";
  dCfg.name = "CylindricalDetector";
  dCfg.builder = containerBuilder;
  dCfg.geoIdGenerator = gig;

  auto detector = DetectorBuilder(dCfg).construct(tContext);

  // Convert the detector
  DetrayConversionUtils::GeometryIdCache geoIdCache;

  detray::io::detector_payload payload =
      DetrayGeometryConverter::convertDetector(geoIdCache, tContext, *detector,
                                               *logger);

  // Test the payload - we have six volumes
  BOOST_CHECK_EQUAL(payload.volumes.size(), 6u);

  // The first volume is the beam pipe volume
  BOOST_CHECK_EQUAL(payload.volumes[0].name, "BeamPipe");
  // The beam pipe should have 1 surface and 5 portals
  // the original cylinder cover is split into 3 portals
  BOOST_CHECK_EQUAL(payload.volumes[0].surfaces.size(), 6u);
  // The second volume should be the negative endcap
  BOOST_CHECK_EQUAL(payload.volumes[1].name, "NegativeEndcap");
  // The negative endcap should have 1 surface and 6 portals
  BOOST_CHECK_EQUAL(payload.volumes[1].surfaces.size(), 7u);
  // Barrel 0,1,2 follow
  BOOST_CHECK_EQUAL(payload.volumes[2].name, "Barrel0");
  BOOST_CHECK_EQUAL(payload.volumes[3].name, "Barrel1");
  BOOST_CHECK_EQUAL(payload.volumes[4].name, "Barrel2");
  // No portal splitting for those, hence we remain at 5 surfaces
  // i.e. 4 portals and one surface
  BOOST_CHECK_EQUAL(payload.volumes[2].surfaces.size(), 5u);
  BOOST_CHECK_EQUAL(payload.volumes[3].surfaces.size(), 5u);
  BOOST_CHECK_EQUAL(payload.volumes[4].surfaces.size(), 5u);
  // Finally the positive endcap
  BOOST_CHECK_EQUAL(payload.volumes[5].name, "PositiveEndcap");
  // The positive endcap should have again 1 surface and 6 portals
  BOOST_CHECK_EQUAL(payload.volumes[5].surfaces.size(), 7u);
}

BOOST_AUTO_TEST_SUITE_END()
