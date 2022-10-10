// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/NavigationState.hpp"
#include "Acts/Experimental/detail/DetectorVolumeLinks.hpp"
#include "Acts/Experimental/detail/NavigationStateUpdators.hpp"
#include "Acts/Experimental/detail/PortalGenerators.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
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

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(CylindricalDetectorVolumePortals) {
  Acts::ActsScalar rInner = 10.;
  Acts::ActsScalar rOuter = 100.;
  Acts::ActsScalar zHalfL = 200.;

  Acts::Transform3 nominal = Acts::Transform3::Identity();

  auto fullCylinderBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(0., rOuter, zHalfL);

  auto portalGenerator = detail::defaultPortalGenerator();

  auto navigationStateUpdator = detail::defaultPortalProvider();

  // Misconfigured - null pointer for bounds
  BOOST_CHECK_THROW(
      DetectorVolumeFactory::construct(
          portalGenerator, tContext, "MisconfiguredFullCylinderVolume", nominal,
          nullptr, detail::defaultPortalProvider()),
      std::invalid_argument);

  // Misconfigured - portal generator not connected
  PortalGenerator unconnected;
  BOOST_CHECK_THROW(
      DetectorVolumeFactory::construct(
          unconnected, tContext, "MisconfiguredFullCylinderVolume", nominal,
          nullptr, detail::defaultPortalProvider()),
      std::invalid_argument);

  // A full cylinder
  auto fullCylinderVolume = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "FullCylinderVolume", nominal,
      std::move(fullCylinderBounds), detail::defaultPortalProvider());

  BOOST_CHECK(fullCylinderVolume ==
              unpackToShared<DetectorVolume>(*fullCylinderVolume));
  BOOST_CHECK(fullCylinderVolume ==
              unpackToShared<const DetectorVolume>(*fullCylinderVolume));

  BOOST_CHECK(fullCylinderVolume->surfaces().size() == 0u);
  BOOST_CHECK(fullCylinderVolume->volumes().size() == 0u);
  BOOST_CHECK(fullCylinderVolume->portals().size() == 3u);

  // A tube cylinder
  auto tubeCylinderBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(rInner, rOuter, zHalfL);

  auto tubeCylinderVolume = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "TubeCylinderVolume", nominal,
      std::move(tubeCylinderBounds), detail::defaultPortalProvider());

  BOOST_CHECK(tubeCylinderVolume->surfaces().size() == 0u);
  BOOST_CHECK(tubeCylinderVolume->volumes().size() == 0u);
  BOOST_CHECK(tubeCylinderVolume->portals().size() == 4u);

  // Let's test the resizing, first inside test: OK
  BOOST_CHECK(tubeCylinderVolume->inside(tContext, Acts::Vector3(50., 0., 0.)));
  // Outside
  BOOST_CHECK(
      not tubeCylinderVolume->inside(tContext, Acts::Vector3(150., 0., 0.)));
  // Resize and check again
  auto biggerBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(0., 200., 200.);
  tubeCylinderVolume->resize(tContext, std::move(biggerBounds),
                             portalGenerator);
  BOOST_CHECK(
      tubeCylinderVolume->inside(tContext, Acts::Vector3(150., 0., 0.)));

  // Check the extent
  auto volumeExtent = tubeCylinderVolume->extent(tContext, 1);
  CHECK_CLOSE_ABS(volumeExtent.min(Acts::binR), 0., 10e-5);
  CHECK_CLOSE_ABS(volumeExtent.max(Acts::binR), 200., 10e-5);
  CHECK_CLOSE_ABS(volumeExtent.min(Acts::binZ), -200., 10e-5);
  CHECK_CLOSE_ABS(volumeExtent.max(Acts::binZ), 200., 10e-5);

  // Misconfigured - wrong bounds
  auto otherBounds =
      std::make_unique<Acts::CuboidVolumeBounds>(100., 100., 100.);
  BOOST_CHECK_THROW(tubeCylinderVolume->resize(tContext, std::move(otherBounds),
                                               portalGenerator),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(CuboidWithCuboid) {
  Acts::ActsScalar bigBox = 100.;
  Acts::ActsScalar smallBox = 10.;

  Acts::Transform3 nominal = Acts::Transform3::Identity();

  auto bigBoxBounds =
      std::make_unique<Acts::CuboidVolumeBounds>(bigBox, bigBox, bigBox);

  auto smallBoxBounds =
      std::make_unique<Acts::CuboidVolumeBounds>(smallBox, smallBox, smallBox);

  auto portalGenerator = detail::defaultPortalGenerator();

  // Create the inner box
  auto innerBox = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "InnerBox", nominal, std::move(smallBoxBounds),
      detail::defaultPortalProvider());

  // A Portal attacher for inner and outer portals
  detail::AllPortalsAttacher ap;
  auto nStateUpdatorStore =
      std::make_shared<detail::NavigationStateUpdator<decltype(ap)>>(
          std::make_tuple(ap));
  NavigationStateUpdator nStateUpdator;
  nStateUpdator.connect<&detail::NavigationStateUpdator<decltype(ap)>::update>(
      nStateUpdatorStore.get());

  ManagedNavigationStateUpdator navStateUpdator;
  navStateUpdator.delegate = std::move(nStateUpdator);
  navStateUpdator.implementation = nStateUpdatorStore;

  // Create the outer box
  std::vector<std::shared_ptr<Acts::Surface>> surfaces = {};
  std::vector<std::shared_ptr<DetectorVolume>> volumes = {innerBox};
  auto outerBox = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "OuterBox", nominal, std::move(bigBoxBounds),
      surfaces, volumes, std::move(navStateUpdator));
  // And connect it to the mother
  detail::setOutsideVolumeLink<DetectorVolume>(*outerBox.get());

  // Check if inside/outside is right: inside big box
  BOOST_CHECK(outerBox->inside(tContext, Acts::Vector3(15., 15., 0.)) == true);
  BOOST_CHECK(outerBox->inside(tContext, Acts::Vector3(5., 5., 0.)) == false);

  // Should sit us just above the internal box
  NavigationState nState;
  outerBox->updateNavigationStatus(nState, tContext, Acts::Vector3(15., 0., 0.),
                                   Acts::Vector3(1., 0., 0), 100., 1.);

  // This should have portals of the outer and inner
  BOOST_TEST(nState.surfaceCandidates.size(), 12u);
  // Current volume is set to outer box
  BOOST_TEST(nState.currentVolume, outerBox.get());

  BOOST_CHECK(nState.surfaceCandidate->objectIntersection.intersection.status ==
              Acts::Intersection3D::Status::reachable);
  CHECK_CLOSE_ABS(
      nState.surfaceCandidate->objectIntersection.intersection.position,
      Acts::Vector3(100, 0, 0), 1e-4);

  // No inside the small box
  nState = NavigationState{};
  outerBox->updateNavigationStatus(nState, tContext, Acts::Vector3(1., 0., 0.),
                                   Acts::Vector3(1., 0., 0), 100., 1.);
  // This has portals only of the inner
  BOOST_TEST(nState.surfaceCandidates.size(), 6u);
  // Current volume is set to inner box
  BOOST_TEST(nState.currentVolume, innerBox.get());
}

BOOST_AUTO_TEST_SUITE_END()
