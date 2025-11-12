// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/interface/IDetectorComponentBuilder.hpp"
#include "Acts/Detector/interface/IGeometryIdGenerator.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"

#include <memory>
#include <stdexcept>
#include <vector>

using namespace Acts;

namespace ActsTests {

class CompBuilder final : public Experimental::IDetectorComponentBuilder {
 public:
  explicit CompBuilder(
      const std::vector<std::shared_ptr<Surface>>& sensitives = {})
      : m_sensitives(sensitives) {}

  Experimental::DetectorComponent construct(
      const GeometryContext& gctx) const final {
    auto bounds = std::make_unique<CuboidVolumeBounds>(10., 10., 10.);
    // Construct the DetectorVolume
    auto dVolume = m_sensitives.empty()
                       ? Experimental::DetectorVolumeFactory::construct(
                             Experimental::defaultPortalGenerator(), gctx,
                             "TestVolume", Transform3::Identity(),
                             std::move(bounds), Experimental::tryAllPortals())
                       : Experimental::DetectorVolumeFactory::construct(
                             Experimental::defaultPortalGenerator(), gctx,
                             "TestVolumeWithSurfaces", Transform3::Identity(),
                             std::move(bounds), m_sensitives, {},
                             Experimental::tryNoVolumes(),
                             Experimental::tryAllPortalsAndSurfaces());

    // Fill the portal container
    Experimental::DetectorComponent::PortalContainer portalContainer;
    for (auto [ip, p] : enumerate(dVolume->portalPtrs())) {
      portalContainer[ip] = p;
    }

    auto geoID = GeometryIdentifier().withVolume(1);
    dVolume->assignGeometryId(geoID);

    return Experimental::DetectorComponent{
        {dVolume},
        portalContainer,
        {{dVolume}, Experimental::tryRootVolumes()}};
  }

 private:
  std::vector<std::shared_ptr<Surface>> m_sensitives;
};

class SurfaceGeoIdGenerator : public Experimental::IGeometryIdGenerator {
 public:
  Experimental::IGeometryIdGenerator::GeoIdCache generateCache() const final {
    return std::any();
  }

  void assignGeometryId(
      Experimental::IGeometryIdGenerator::GeoIdCache& /*cache*/,
      Experimental::DetectorVolume& dVolume) const final {
    for (auto [is, s] : enumerate(dVolume.surfacePtrs())) {
      auto geoID = GeometryIdentifier().withSensitive(is + 1);
      s->assignGeometryId(geoID);
    }
  }

  void assignGeometryId(
      Experimental::IGeometryIdGenerator::GeoIdCache& /*cache*/,
      Experimental::Portal& /*portal*/) const final {}

  void assignGeometryId(
      Experimental::IGeometryIdGenerator::GeoIdCache& /*cache*/,
      Surface& /*surface*/) const final {}
};

GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(DetectorSuite)

BOOST_AUTO_TEST_CASE(DetectorBuilder_Misconfigured) {
  // Detector builder
  Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxiliary = "*** Test X * Misconfigued ***";
  dCfg.name = "EmptyCylinder";
  dCfg.builder = nullptr;

  BOOST_CHECK_THROW(auto a = Experimental::DetectorBuilder(dCfg),
                    std::invalid_argument);

  // Detector builder with sensitives but no assigned geometry ids to them
  ActsTests::DetectorElementStub detElement0(
      Transform3::Identity(), std::make_shared<RectangleBounds>(5., 5.), 0.1);
  ActsTests::DetectorElementStub detElement1(
      Transform3::Identity(), std::make_shared<RectangleBounds>(5., 5.), 0.1);

  std::vector<std::shared_ptr<Surface>> sensitives;
  sensitives.push_back(detElement0.surface().getSharedPtr());
  sensitives.push_back(detElement1.surface().getSharedPtr());
  dCfg.builder = std::make_shared<CompBuilder>(sensitives);

  BOOST_CHECK_THROW(Experimental::DetectorBuilder(dCfg).construct(tContext),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(DetectorBuilder_test) {
  // Detector builder
  Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxiliary = "*** Test : Detector ***";
  dCfg.name = "TestDetector";
  dCfg.builder = std::make_shared<CompBuilder>();

  BOOST_CHECK_NO_THROW(auto a = Experimental::DetectorBuilder(dCfg));

  auto detector = Experimental::DetectorBuilder(dCfg).construct(tContext);

  BOOST_CHECK_EQUAL(detector->name(), "TestDetector");
  BOOST_CHECK_EQUAL(detector->rootVolumes().size(), 1);
}

BOOST_AUTO_TEST_CASE(DetectorBuilder_testWithSurfaces) {
  // Detector builder
  Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxiliary = "*** Test : Detector ***";
  dCfg.name = "TestDetectorWithSurfaces";
  dCfg.builder = std::make_shared<CompBuilder>();

  // Test detector with surfaces
  ActsTests::DetectorElementStub detElement0(
      Transform3::Identity(), std::make_shared<RectangleBounds>(5., 5.), 0.1);
  ActsTests::DetectorElementStub detElement1(
      Transform3::Identity(), std::make_shared<RectangleBounds>(5., 5.), 0.1);

  std::vector<std::shared_ptr<Surface>> sensitives;
  sensitives.push_back(detElement0.surface().getSharedPtr());
  sensitives.push_back(detElement1.surface().getSharedPtr());
  dCfg.builder = std::make_shared<CompBuilder>(sensitives);
  dCfg.geoIdGenerator = std::make_shared<SurfaceGeoIdGenerator>();

  auto detector = Experimental::DetectorBuilder(dCfg).construct(tContext);
  BOOST_CHECK_EQUAL(detector->name(), "TestDetectorWithSurfaces");
  BOOST_CHECK_EQUAL(
      detector->volumes()[0]->surfaces()[0]->geometryId().sensitive(), 1);
  BOOST_CHECK_EQUAL(
      detector->volumes()[0]->surfaces()[1]->geometryId().sensitive(), 2);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
