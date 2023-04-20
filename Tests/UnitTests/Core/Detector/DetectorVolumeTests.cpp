// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/PortalHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/DetectorVolumeUpdators.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <exception>
#include <memory>

/// Unpack to shared - simply to test the getSharedPtr mechanism
///
/// @tparam referenced_type is the type of the referenced object
///
/// @param rt is the referenced object
///
/// @returns a shared pointer
template <typename referenced_type>
std::shared_ptr<referenced_type> unpackToShared(referenced_type& rt) {
  return rt.getSharedPtr();
}

using namespace Acts::Experimental;

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(CylindricalDetectorVolumePortals) {
  Acts::ActsScalar rInner = 10.;
  Acts::ActsScalar rOuter = 100.;
  Acts::ActsScalar zHalfL = 200.;

  Acts::Transform3 nominal = Acts::Transform3::Identity();

  auto fullCylinderBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(0., rOuter, zHalfL);

  auto portalGenerator = defaultPortalGenerator();

  // Misconfigured - null pointer for bounds
  BOOST_CHECK_THROW(
      DetectorVolumeFactory::construct(portalGenerator, tContext,
                                       "MisconfiguredFullCylinderVolume",
                                       nominal, nullptr, tryAllPortals()),
      std::invalid_argument);

  // Misconfigured - portal generator not connected
  PortalGenerator unconnected;
  BOOST_CHECK_THROW(
      DetectorVolumeFactory::construct(unconnected, tContext,
                                       "MisconfiguredFullCylinderVolume",
                                       nominal, nullptr, tryAllPortals()),
      std::invalid_argument);

  // A full cylinder
  auto fullCylinderVolume = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "FullCylinderVolume", nominal,
      std::move(fullCylinderBounds), tryAllPortals());

  BOOST_CHECK(fullCylinderVolume ==
              unpackToShared<DetectorVolume>(*fullCylinderVolume));
  BOOST_CHECK(fullCylinderVolume ==
              unpackToShared<const DetectorVolume>(*fullCylinderVolume));

  BOOST_CHECK(fullCylinderVolume->surfaces().empty());
  BOOST_CHECK(fullCylinderVolume->volumes().empty());
  BOOST_CHECK(fullCylinderVolume->portals().size() == 3u);

  // A tube cylinder
  auto tubeCylinderBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(rInner, rOuter, zHalfL);

  auto tubeCylinderVolume = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "TubeCylinderVolume", nominal,
      std::move(tubeCylinderBounds), tryAllPortals());

  BOOST_CHECK(tubeCylinderVolume->surfaces().empty());
  BOOST_CHECK(tubeCylinderVolume->volumes().empty());
  BOOST_CHECK(tubeCylinderVolume->portals().size() == 4u);

  // Let's test the resizing, first inside test: OK
  BOOST_CHECK(tubeCylinderVolume->inside(tContext, Acts::Vector3(50., 0., 0.)));
  // Outside
  BOOST_CHECK(
      not tubeCylinderVolume->inside(tContext, Acts::Vector3(150., 0., 0.)));

  // Check the extent
  auto volumeExtent = tubeCylinderVolume->extent(tContext, 1);
  CHECK_CLOSE_ABS(volumeExtent.min(Acts::binR), 10., 10e-5);
  CHECK_CLOSE_ABS(volumeExtent.max(Acts::binR), 100., 10e-5);
  CHECK_CLOSE_ABS(volumeExtent.min(Acts::binZ), -200., 10e-5);
  CHECK_CLOSE_ABS(volumeExtent.max(Acts::binZ), 200., 10e-5);
}

BOOST_AUTO_TEST_CASE(UpdatePortal) {
  Acts::Transform3 nominal = Acts::Transform3::Identity();

  auto fullCylinderBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(0., 10., 100.);

  auto portalGenerator = defaultPortalGenerator();

  auto fullCylinderVolume = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "FullCylinderVolume", nominal,
      std::move(fullCylinderBounds), tryAllPortals());

  auto cylinderSurface =
      Acts::Surface::makeShared<Acts::CylinderSurface>(nominal, 10., 100.);

  auto cylinderPortal = Acts::Experimental::Portal::makeShared(cylinderSurface);

  fullCylinderVolume->updatePortal(cylinderPortal, 2u);

  BOOST_CHECK(fullCylinderVolume->portals()[2u] == cylinderPortal.get());
}

BOOST_AUTO_TEST_CASE(CuboidWithCuboid) {
  Acts::ActsScalar bigBox = 100.;
  Acts::ActsScalar smallBox = 10.;

  Acts::Transform3 nominal = Acts::Transform3::Identity();

  auto bigBoxBounds =
      std::make_unique<Acts::CuboidVolumeBounds>(bigBox, bigBox, bigBox);

  auto smallBoxBounds =
      std::make_unique<Acts::CuboidVolumeBounds>(smallBox, smallBox, smallBox);

  auto portals = defaultPortalGenerator();
  auto generatePortalsUpdateInternals = defaultPortalAndSubPortalGenerator();

  // Create the inner box
  auto innerBox = DetectorVolumeFactory::construct(
      portals, tContext, "InnerBox", nominal, std::move(smallBoxBounds),
      tryAllPortals());

  std::vector<std::shared_ptr<Acts::Surface>> surfaces = {};
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> volumes = {
      innerBox};

  // Create the outer box and insert the inner box, use a portal generator
  // with sub portal registration
  auto outerBox = DetectorVolumeFactory::construct(
      generatePortalsUpdateInternals, tContext, "OuterBox", nominal,
      std::move(bigBoxBounds), surfaces, volumes, tryAllSubVolumes(),
      tryAllPortals());

  // Check that we are within the outer box
  Acts::Experimental::NavigationState nState;
  nState.position = Acts::Vector3(-50., 5., 0.);
  nState.direction = Acts::Vector3(1., 0., 0.);

  BOOST_CHECK(outerBox->inside(tContext, nState.position));
  nState.currentVolume = outerBox.get();

  outerBox->updateNavigationState(tContext, nState);

  // We should have 12 candidates, 6 inner, 6 outer portals
  BOOST_CHECK(nState.surfaceCandidates.size() == 12u);
}

BOOST_AUTO_TEST_SUITE_END()
